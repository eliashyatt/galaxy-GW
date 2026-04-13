#!/bin/bash
export IFS=','

fin='O3a_events_HLV_SEOBNRv4PHM_maxL.csv'
#awk -F, 'NR==4{print $2}' O3a_events_HLV_C01:SEOBNRv4PHM_maxL.csv
# cat O3a_events_HLV_SEOBNRv4PHM_maxL.csv | awk -F, 'NR==4{print $2}'

#event_index=3
# loop the events
event_index=2
for((event_index=2;event_index<=20;event_index++))
do
# extracting the MaxLparameters used in lalapps
event=`cat $fin | awk -F, 'NR==1{print $'$event_index'}'`
geocent_time=`cat $fin | awk -F, 'NR==2{printf "%0.3f\n", $'$event_index'}'`
luminosity_distance=`cat $fin | awk -F, 'NR==4{printf "%0.3f\n", $'$event_index'}'`
m1=`cat $fin | awk -F, 'NR==5{printf "%0.3f\n", $'$event_index'}'`
m2=`cat $fin | awk -F, 'NR==6{printf "%0.3f\n", $'$event_index'}'`
dec=`cat $fin | awk -F, 'NR==7{printf "%0.3f\n", $'$event_index'}'`
ra=`cat $fin | awk -F, 'NR==8{printf "%0.3f\n", $'$event_index'}'`
theta_jn=`cat $fin | awk -F, 'NR==9{printf "%0.3f\n", $'$event_index'}'`
psi=`cat $fin | awk -F, 'NR==11{printf "%0.3f\n", $'$event_index'}'`
phase=`cat $fin | awk -F, 'NR==12{printf "%0.3f\n", $'$event_index'}'`
a_1=`cat $fin | awk -F, 'NR==13{printf "%0.3f\n", $'$event_index'}'`
a_2=`cat $fin | awk -F, 'NR==14{printf "%0.3f\n", $'$event_index'}'`

echo $event,$geocent_time
 

# running lalapp to make injections
# covert from radian to degree 1rad=57.3deg

time=${geocent_time%.*}
gps_start=$(echo "$time-210"|bc)
#gps_start=$(echo "$time-190"|bc)
gps_end=$(echo "$time+50"|bc)
echo $gps_start,$time,$gps_end

theta_jn_deg=$(echo "scale=3;$theta_jn*57.3/1"|bc)
phase_deg=$(echo "scale=3;$phase*57.3/1"|bc)
psi_deg=$(echo "scale=3;$psi*57.3/1"|bc)

spin1_max=$(printf "%.3f" `echo "scale=4;($a_1+0.001)/1"|bc`)
spin2_max=$(printf "%.3f" `echo "scale=4;($a_2+0.001)/1"|bc`)

luminosity_distance=$(echo "($luminosity_distance*1000)/1"|bc)
luminosity_distance_max=$(echo "($luminosity_distance+1)/1"|bc)


echo $theta_jn_deg,$phase_deg
echo $psi_deg,$gps_start,$gps_end
echo $spin1_max,$spin2_max,$luminosity_distance_max

# check if ra_deg is in range[-180,180] as lalapp requirs
ra_deg_tmp=$(echo "scale=3;$ra*57.3/1"|bc)
ra_deg_tmp_int=${ra_deg_tmp%.*}

if [ $ra_deg_tmp_int -gt 180 ]
then
    	echo "ra is greate than 180"
        ra_deg=$(echo "$ra_deg_tmp-360"|bc)
else
    	ra_deg=$ra_deg_tmp
fi
echo $ra_deg_tmp,$ra_deg


# check	if dec_deg is in range[-90,90]as lalapp requirs
dec_deg_tmp=$(echo "scale=3;$dec*57.3/1"|bc)
dec_deg_tmp_int=${dec_deg_tmp%.*}

if [ $dec_deg_tmp_int -gt 90 ]
then
    	echo "dec is greate than 90"
        dec_deg=$(echo "$dec_deg_tmp-360"|bc)
else
        dec_deg=$dec_deg_tmp
fi
echo $dec_deg_tmp,$dec_deg


## run lalapp to get the xml file

echo lalapps_inspinj --gps-start-time $gps_start --gps-end-time $gps_end --time-step 200 \
--m-distr fixMasses --fixed-mass1 $m1 --fixed-mass2 $m2 \
--i-distr fixed --fixed-inc $theta_jn_deg --polarization $psi_deg --fixed-coa-phase $phase_deg \
--l-distr fixed --longitude $ra_deg  --latitude $dec_deg \
--waveform SEOBNRv4PHM --f-lower 10 --taper-injection startend \
--d-distr uniform --min-distance $luminosity_distance --max-distance $luminosity_distance_max \
--enable-spin --min-spin1 $a_1 --max-spin1 $spin1_max --min-spin2 $a_2 --max-spin2 $spin2_max \
--band-pass-injection

lalapps_inspinj --gps-start-time $gps_start --gps-end-time $gps_end --time-step 200 \
--m-distr fixMasses --fixed-mass1 $m1 --fixed-mass2 $m2 \
--i-distr fixed --fixed-inc $theta_jn_deg --polarization $psi_deg --fixed-coa-phase phase_deg \
--l-distr fixed --longitude $ra_deg  --latitude $dec_deg \
--waveform SEOBNRv4PHM --f-lower 10 --taper-injection startend \
--d-distr uniform --min-distance $luminosity_distance --max-distance $luminosity_distance_max \
--enable-spin --min-spin1 $a_1 --max-spin1 $spin1_max --min-spin2 $a_2 --max-spin2 $spin2_max \
--band-pass-injection 
echo "lalapp simulation is over"
echo "start check the xml file"
ls *xml
#lal_output='HL-INJECTIONS_1-'$time'-0.xml'
lal_output=`ls HL-INJECTIONS_1*`
echo $lal_output
lal_new_name='./output_xml_sh/'$event'_10s_before_lalapp.xml'
echo $lal_new_name

mv $lal_output $lal_new_name

coinc_name='./output_xml_sh/bayestar_coinc_'$event'_10s_before_lalapp.xml'

bayestar-realize-coincs -o $coinc_name $lal_new_name --reference-psd ./psd/psd_O3_actual_LHV.xml --detector H1 L1 V1 \
--measurement-error gaussian-noise --snr-threshold 4.0 --net-snr-threshold 8.0 --min-triggers 2 --keep-subthreshold

export OMP_NUM_THREADS=4
bayestar-localize-coincs $coinc_name --output ./sm_tmp --waveform SEOBNRv4PHM --f-low 10


fits_name='./skymaps/'$event'_10s_before_lalapp.fits'
mv ./sm_tmp/1.fits $fits_name
done

python read_all_skymaps.py