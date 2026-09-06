import pandas as pd
import subprocess
import os
import numpy as np
import argparse

# Input CSV file
#fin = 'inj_parameters_15000.csv'    #input

def generate_skymaps(fin, max_count, foldername, psdname):

    # Create output directories if they don't exist
    #os.makedirs('./output_xml_sh', exist_ok=True)
    #os.makedirs('./skymaps', exist_ok=True)
    #os.makedirs('./sm_tmp', exist_ok=True)

    # Read CSV into a DataFrame
    df = pd.read_csv(fin)
    print(df.columns)

    #if df.columns[0].startswith("Unnamed"):
        #df.rename(columns={df.columns[0]: "event"}, inplace=True)

    df_fail = pd.DataFrame(columns=df.columns)
    failed_rows = []
    count = 0    #use if continuing from a previous run, otherwise set to 0; this is the number of skymaps already generated
    num = 0
    #max_count = 3000

    # Loop over events (rows)
    for _, row in df.iterrows():
        num += 1
        if num < 0:    #use if continuing from a previous run, otherwise set to 0; this is the event number of the last skymap generated
            continue
    
        # Extract variables
        event = f"{int(row['event']):05d}"
        geocent_time = float(row['geocent_time']) #distribution
        luminosity_distance = float(row['luminosity_distance']) #from halos
        m1 = float(row['mass_1_source']) #distribution
        m2 = float(row['mass_2_source']) #distribution from mass ratio q=m2/m1
        dec = float(row['dec']) #from halos
        ra = float(row['ra']) #from halos
        theta_jn = float(row['theta_jn']) #random
        psi = float(row['psi']) #random
        phase = float(row['phase']) #random
        a_1 = float(row['a_1']) #distribution
        a_2 = float(row['a_2']) #distribution

        print(f"{event},{geocent_time:.3f}")

        # GPS times
        time = int(geocent_time)
        gps_start = time
        gps_end = time + 10

        # Convert radians to degrees
        theta_jn_deg = theta_jn * 57.3
        phase_deg = phase * 57.3
        psi_deg = psi * 57.3

        # Spin
        spin1_max = round(a_1 + 0.001, 3)
        spin2_max = round(a_2 + 0.001, 3)

        # Distances
        luminosity_distance_km = int(luminosity_distance * 1000)
        luminosity_distance_max = luminosity_distance_km + 1

        # RA range check / conversion
        ra_deg = ra * 180 / np.pi
        if ra_deg > 180:
            ra_deg -= 360

        # DEC range check / conversion
        dec_deg = dec * 180 / np.pi
        #if int(dec_deg) >= 90:
            #dec_deg -= 360

        # Run lalapps_inspinj
        lal_cmd = [
            "lalapps_inspinj",
            "--gps-start-time", str(gps_start),
            "--gps-end-time", str(gps_end),
            "--time-step", "20",
            #"--time-interval", "50",
            "--m-distr", "fixMasses",
            "--fixed-mass1", str(m1),
            "--fixed-mass2", str(m2),
            "--i-distr", "fixed",
            "--fixed-inc", str(theta_jn_deg),
            "--polarization", str(psi_deg),
            "--fixed-coa-phase", str(phase_deg),
            "--l-distr", "fixed",
            "--longitude", str(ra_deg),
            "--latitude", str(dec_deg),
            "--waveform", "IMRPhenomXPHM-SpinTaylor",
            "--f-lower", "10",
            "--taper-injection", "startend",
            "--d-distr", "uniform",
            "--min-distance", str(luminosity_distance_km),
            "--max-distance", str(luminosity_distance_max),
            "--enable-spin",
            "--min-spin1", str(a_1),
            "--max-spin1", str(spin1_max),
            "--min-spin2", str(a_2),
            "--max-spin2", str(spin2_max),
            "--band-pass-injection"
        ]
        print(f"Running lalapps_inspinj for {event}...")
        subprocess.run(lal_cmd, check=True)

        # Rename XML output
        lal_output_files = [f for f in os.listdir('.') if f.startswith("HL-INJECTIONS_1") and f.endswith(".xml")]
        if not lal_output_files:
            raise FileNotFoundError(f"No lalapps XML output found for {event}")
        lal_output = lal_output_files[0]
        lal_new_name = f'./output_xml_sh/{foldername}/{event}.xml'
        os.rename(lal_output, lal_new_name)

        # Run bayestar
        coinc_name = f'./output_xml_sh/{foldername}/bayestar_coinc_{event}.xml'
        bayestar_realize_cmd = [
            "bayestar-realize-coincs",
            "-o", coinc_name,
            lal_new_name,
            "--reference-psd", f"./psd/{psdname}.xml",
            #"--detector", "H1", "L1", "V1", #L1:aLIGOaLIGODesignSensitivityT1800044 #H1:aLIGOaLIGODesignSensitivityT1800044 #V1:AdVDesignSensitivityP1200087
            "--detector", "H1", "L1", "V1",
            "--measurement-error", "gaussian-noise",
            "--snr-threshold", "4.0",
            "--net-snr-threshold", "8.0",
            "--min-triggers", "3",
            "--keep-subthreshold"
        ]
        try:
            subprocess.run(bayestar_realize_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Bayestar failed for {event} with exit code {e.returncode}")
            print("stdout:", e.stdout)
            print("stderr:", e.stderr)
            continue  # skip to next event

        os.environ["OMP_NUM_THREADS"] = "4"
        bayestar_localize_cmd = [
            "bayestar-localize-coincs",
            coinc_name,
            "--output", "./sm_tmp",
            "--waveform", "IMRPhenomXPHM-SpinTaylor",
            "--f-low", "10",
        ]
        subprocess.run(bayestar_localize_cmd, check=True)

        # Move final FITS file
        fits_name = f'./skymaps/{foldername}/{event}_{psdname}.fits'
        fits_files = [f for f in os.listdir('./sm_tmp') if f.endswith('.fits')]

        if fits_files:
            fits_input = os.path.join('./sm_tmp', fits_files[0])
            os.rename(fits_input, fits_name)
            print(f"Saved skymap: {fits_name}")
            count = count + 1
        else:
            print(f"No FITS file found for {event}")
            failed_rows.append(row)

        if count == max_count:
            break
        
    df_fail = pd.DataFrame(failed_rows)
    #df_fail.to_csv("failed_rows_1000_plus.csv", index=False)
    df_fail.to_csv(f'failed_rows.csv', index=False)

    success_rows_df = df[~df["event"].isin(df_fail["event"])]
    success_rows_df = success_rows_df[:max_count]

    #success_rows_df.to_csv("success_rows_1000_plus.csv", index=False)
    success_rows_df.to_csv(f'success_rows_{max_count}.csv', index=False)

    # Run final Python script
    #subprocess.run(["python", "read_all_skymaps.py"], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates skymaps based on injection parameters")
    parser.add_argument("--input_file", help="CSV file of injection parameters", required=True)
    parser.add_argument("--n_skymaps", help="Number of wanted skymaps to be generated", required=True)
    parser.add_argument("--foldername", help="Name of folder inside /skymaps/ to save skymaps to", required=True)
    parser.add_argument("--psdname", help="Name of the psd inside /psd/ being used WITHOUT .xml", required=True)

    args = parser.parse_args()

    fin = args.input_file
    max_count = int(args.n_skymaps)
    foldername = args.foldername
    psdname = args.psdname

    os.makedirs(f'./skymaps/{foldername}', exist_ok=True)

    os.makedirs(f'./output_xml_sh/{foldername}', exist_ok=True)
    #os.makedirs('./skymaps', exist_ok=True)
    os.makedirs(f'./sm_tmp/', exist_ok=True)

    generate_skymaps(fin, max_count, foldername, psdname)

    

    