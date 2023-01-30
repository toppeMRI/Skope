# Pulse sequences for Skope field camera installation and use on GE scanners

For general scan instructions with TOPPE, 
see README.md in <https://github.com/jfnielsen/TOPPEpsdSourceCode/>.


## Setup

In the following, the working directory is assumed to be ~/Skope/GE/.

Install this toolbox, and the TOPPE and MIRT toolboxes:
```
$ pwd
~/Skope/GE/
$ git clone git@github.com:toppeMRI/Skope.git
$ git clone git@github.com:toppeMRI/toppe.git
$ git clone git@github.com:JeffFessler/mirt.git (optional -- for image display)
$ cd Skope
```

In Matlab:
```
>> pwd
~/Skope/GE/Skope/
>> addpath ../toppe/
>> cd ../mirt; setup; cd ../Skope;
```

## Run the 3D FLASH sequence to verify correct installation of the TOPPE interpreter

1. Create scan files (flash.tar):
```
>> pwd
~/Skope/GE/Skope/
>> write3dflash;  
```

2. Copy the contents of flash.tar to /usr/g/research/pulseq/flash/
on the scanner host computer.

3. Move toppe0.entry to /usr/g/research/pulseq/:
```
$ pwd
/usr/g/research/pulseq/flash/
$ mv toppe0.entry ..
```

4. Prescribe the TOPPE interpreter with the following settings:
```
opuser1 = 0   (selects toppe0.entry as the entry point)
oprbw = 31.25 
```
Other settings are as shown in the 
[EPIC source code repository](https://github.com/jfnielsen/TOPPEpsdSourceCode/).

5. (optional) Reconstruct and display:
```
>> recon3dflash;
```


## Run the Skope field camera calibration sequence

Enter Manual Prescan, then click 'scan TR'.

This sequence outputs TTL pulses on the J6 BNC out in the system rack.


## Run the spiral sequence and measure the gradients with the Skope field camera

1. Create the spiral scan files (spiral.tar):
```
>> writespiral;
```

2. Copy the contents of spiral.tar to /usr/g/research/pulseq/spiral/ 
on the scanner host computer.

3. Move toppe1.entry to /usr/g/research/pulseq/:
```
$ pwd
/usr/g/research/pulseq/spiral/
$ mv toppe1.entry ..
```

4. Select the spiral scan:
   Select the `Research/display CVs' menu, and set opuser1 = 1. 
   This will select toppe1.entry as the new entry point for the interpreter.

5. Select `Download`, and click `Scan`.

6. This sequence outputs 1ms TTL pulses at regular intervals.
You may edit writespiral.m and modules.txt to control the location
and frequency of the TTL pulses.

7. For additional details on using the Skope field camera, see
<https://github.com/SkopeMagneticResonanceTechnologies/Pulseq-Sequences>.
