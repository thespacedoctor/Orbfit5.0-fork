# Orbfit v5.0.5 Fork

This is a 'fork' of the [Orbfit code found here](http://adams.dm.unipi.it/orbfit/). At the time I forked the code it was at v5.0.5.

The main contribution I've made is to craft an new executable called `ephem` that will generate an ephemeris for a single object, at a single epoch given a specific observatory code. This hacked/new code can be found in `/src/orbfit/ephem.f90`

I've also made a few minor modifications to the main body of code, including to increase the precision of many of the ephemeris outputs.

## Installation

To install this fork of Orbfit, you should be able to clone this repo to your machine and then [follow these install instructions](http://psweb.mp.qub.ac.uk/dry//blog/2017/09/15/Installing-OrbFit-5.0-on-macOS.html) but using the code in this repo and not the code downloaded from the Orbfit distribution site.

```bash
git clone git@github.com:thespacedoctor/Orbfit5.0-fork.git
```

## ephem Usage

```bash
usage:
   ephem <obscode> <mjd> <objectName>

   <obscode>:        observatory code (use 500 for geocentric)
   <mjd>:            the modified julian date of the ephemeris you wish to generate (UTC)
   <objectName>:     the ID of the asteroid you wish to generate an ephemeris for (MPC number or name)

cmdline options:

  -h, --help        print usage information and exit
```

For example, to generate an ephemeris for asteroid `10547` from the ATLAS observatory at Haleakela on MJD `57916.0` run:

```bash
$ ephem T05 57916.0 10547
mjd,ra_deg,dec_deg,apparent_mag,observer_distance,heliocentric_distance,sun_obs_target_angle,phase_angle,galactic_latitude,ra_arcsec_per_hour,dec_arcsec_per_hour,object_name,obscode
57916.000000,15.077432,2.039571,18.123,2.1104503,1.9411821,66.3842,28.6332,-60.7539,65.8674,26.0504,10547,T05
```