# Orbfit v5.0.5 Fork

This is a 'fork' of the [Orbfit code found here](http://adams.dm.unipi.it/orbfit/). At the time I forked the code it was at v5.0.5.

The main contribution I've made is to craft an new executable called `ephem` that will generate an ephemeris for a single object, at a single epoch given a specific observatory code. This hacked/new code can be found in `/src/orbfit/ephem.f90`

I've also made a few minor modifications to the main body of code, including to increase the precision of many of the ephemeris outputs.

## Installation

To install this fork of Orbfit, you should be able to clone this repo to your machine and then [follow these install instructions](http://psweb.mp.qub.ac.uk/dry//blog/2017/09/15/Installing-OrbFit-5.0-on-macOS.html) but using the code in this repo and not the code downloaded from the Orbfit distribution site.

```bash
git clone git@github.com:thespacedoctor/Orbfit5.0-fork.git
```

