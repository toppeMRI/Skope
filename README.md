# Pulse sequences for assisting with Skope field camera calibration on GE scanners

## 3D FLASH sequence for testing TOPPE interpreter

### Data acquisition

Create scan files:
```
>> write3dflash; 
```

To scan, follow the instructions in <https://github.com/jfnielsen/TOPPEpsdSourceCode/>.

Additional scan instructions:  
* oprbw = 31.25

### Reconstruct:
```
>> recon3dflash;
```

