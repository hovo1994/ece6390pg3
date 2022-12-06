# satellitetracking.py 

you can stuff the dependencies into a requirements.txt and insatll it or do it manually via pip in a virtual environment 

## Create Virtual Env

```
python3 -m venv myvenv
source myvenv/bin/activate
pip install -r requirments.txt
```
OR manually
```
pip install astropy poliastro matplotlib geopy
```
## Dependencies

Note this installation might take a little while, but once installed you can run the ground track simulator. 

```
astropy
poliastro
matplotlib
geopy
```



# plotdata.py

follow same venv instructions from above if you just want to plot data, then only install matplotlib.. 

## Dependencies 

`pip install matplotlib`

```
(satcomms) hov@FancyComputingBox FinalPaper % python plotdata.py --help
usage: Plot stuff [-h] [-f FILENAME]

plots csv stuff

optional arguments:
  -h, --help            show this help message and exit
  -f FILENAME, --filename FILENAME
                        csv file name
```

plotdata.py -f ground_track_v_10_with_rcv_pwr.csv