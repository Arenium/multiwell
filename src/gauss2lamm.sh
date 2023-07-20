#!/bin/bash
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                              GAUSS2LAMM-2010.0
#
#                                 Jan 2010
#
#                             T. J. Dhilip Kumar
#                             thogluva@umich.edu
#
#                           University of Michigan
#                       Ann Arbor, Michigan 48109-2143
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo Enter Gaussian Output Filename
read fname

# Checking keyword NOSYM||Z-MATRIX in the output
dos2unix -q  $fname
dd if=$fname of=f1name conv=ucase 
a="$( grep -o -m 1 " NOSYM" f1name)"  
b=" NOSYM"
aa="$( grep -o -m 1 " Z-MATRIX" f1name)"  
bb=" Z-MATRIX"
rm -f f1name 

if [[ "$b" == "$a" || "$bb" == "$aa" ]]; then
   if [ "$b" == "$a" ]; then
     echo
     echo Keyword $a found in the output
   else
     if [ "$bb" == "$aa" ]; then
       echo 
       echo Keyword $aa found in the output
     fi 
   fi
   echo

#Calculation information
echo Enter total number of atoms in the molecule 
read TNA
echo
echo Enter MIN angle,  No. of Points  and Stepsize in degrees used in the calculation
read AMIN NOP AINT

#Calculating Maximum angle
AMAX=$(echo $AINT*$NOP | bc)

#Printing the calculation information
echo "[Title]                                       "  > lamm.dat   
echo "Calculated at [      ] level of theory      "  >> lamm.dat
echo                                                 >> lamm.dat
echo "$TNA               !Number of atoms in the molecule"  >> lamm.dat         
echo                                                 >> lamm.dat
echo "$AMIN, $NOP, $AINT       ! Minumum angle, No. of Points and Stepsize (in degrees)" >> lamm.dat
echo                                                 >> lamm.dat
#echo "Replace atomic no. with atomic masses (amu) for the elements: "  >> lamm.dat

#Variable initialization
UNIT=1
UNIT4=4
TNA1=$((TNA + $UNIT4))
TNA2=$((TNA + $UNIT))
TNA3=$((TNA - $UNIT))

# Collecting orientation data from the output
   c="$( grep -o -m 1 "Input orient" $fname)"
   d="$( grep -o -m 1 "Z-Matrix orien" $fname)"
   e="Input orient"
   f="Z-Matrix orien"
   if [ "$e" == "$c" ]; then
      echo
      echo "Proceeding to collect data from $fname"
      awk '/Input ori/{s=x}{s=s$0"\n"}/Optimized Pa/{print s}' $fname > t1
      grep -A $TNA1  "Input orient" t1 > t2
      grep -A $TNA2  " X           Y           Z" t2 > t3
      grep -B $TNA3  "   $TNA   " t3 > t4
      grep -m 1 -B $TNA3 " $TNA " t3 | awk '{print $2}' >> lamm.dat
      echo >> lamm.dat
      awk '{print $4, $5, $6}' t4 >> lamm.dat   
      echo >> lamm.dat

# Collecting energy data from the output
      grep "SCF Done:" t1 | awk '{print $5}'  > lamme.dat
      rm -f t1 t2 t3 t4 

   else
      if [ "$f" == "$d" ]; then
      echo
      echo "Proceeding to collect data from $fname"
      awk '/Z-Matrix orient/{s=x}{s=s$0"\n"}/Optimized Pa/{print s}' $fname > t1
      grep -A $TNA1 "Z-Matrix orient" t1 > t2
      grep -A $TNA2 " X           Y           Z" t2 > t3
      grep -B $TNA3 "   $TNA   " t3 > t4
      grep -m 1 -B $TNA3 " $TNA " t3 | awk '{print $2}' >> lamm.dat
      echo >> lamm.dat
      awk '{print $4, $5, $6}' t4 >> lamm.dat   
      echo >> lamm.dat

# Collecting energy data from the output
      grep "SCF Done:" t1 | awk '{print $5}'  > lamme.dat
      rm -f t1 t2 t3 t4  
      fi
   fi

   rm -f lamm_erg.dat lamm_ergcm.dat
   i=$AMIN
   while [ $i -le $AMAX ]
     do
       echo $i >> lammang.dat
       i=`expr $i + $AINT`
     done

# Converting energy into wavenumbers with arbitrary zero minimum
   cp lamme.dat sort.dat
   s="$(sort sort.dat | tail -1)"
   echo
   echo 'Lowest energy in the surface scan:' $s 'hartrees'
   rm -f sort.dat

   echo "'WARNING: The zero of the relative energy is arbitrary minimum'" > lamm_ergcm.dat
   echo >> lamm_ergcm.dat
 
   cm=219474.63
   rm -f ecm.dat
   while read ehart 
   do
   ediff=$(echo $ehart - $s | bc)
   echo $ediff*$cm | bc >> ecm.dat
   done < "lamme.dat"

   awk '{ getline ln < "lammang.dat"; print ln" "$0 }' ecm.dat >> lamm_ergcm.dat
   cat lamm.dat lamm_ergcm.dat  > lamm2.dat
   rm -f lammang.dat lamme.dat ecm.dat lamm_ergcm.dat
   rm -f lamm.dat 
   mv lamm2.dat lamm.dat

else
echo "Keyword NOSYM or Z-MATRIX not found in the output"
echo "Scan using the keyword NOSYM or Z-MATRIX in GAUSSIAN" 
fi

