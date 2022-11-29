
from astropy import time, units as u
import argparse
from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from math import pi, sqrt, tan, acos, atan2, cos, sin
from poliastro.czml.extract_czml import CZMLExtractor
from poliastro.earth.plotting import GroundtrackPlotter
import sys
import geopy.distance

# Earth focused modules, ISS example orbit and time span generator
from poliastro.earth import EarthSatellite
from poliastro.earth.plotting import GroundtrackPlotter
from poliastro.util import time_range
import matplotlib.pyplot   as     plt

RADIUS_EARTH = 6380000 # KM

# NOTES: Use 

class SphericalCoordinates:
    """Position in spherical coordinates"""
    def __init__(self, r, theta, phi, x, y, z):
        self.r = r
        self.theta = theta
        self.phi = phi
        self.x = x
        self.y = y
        self.z = z
        # self.ssp_x = None
        # self.ssp_y = None
        # self.ssp_z = None
        self.altitude = self.r - RADIUS_EARTH

        # self._computeSSPCartesian()

    @classmethod
    def from_cartesian(cls, x, y, z):
        x = x * 1000
        y = y * 1000
        z = z * 1000
        r       =  sqrt(x*x + y*y + z*z)
        theta   =  acos(z/r)*180/ pi
        phi     =  atan2(y,x)*180/ pi
        return cls(r, theta, phi, x, y, z)

    # def _computeSSPCartesian(self):
    #     r_ssp = RADIUS_EARTH
    #     self.ssp_x = r_ssp * sin(self.phi)*cos(self.theta)
    #     self.ssp_y = r_ssp * sin(self.phi)*sin(self.theta)
    #     self.ssp_z = r_ssp * cos(self.phi)        





class SatPosEpoch:
    """Compute Satellite Beam from position and frequency"""
    def __init__(self, Position, lat, lon, freq=5.8e9, D_T=333):
        self.Position = Position
        self.D_T = D_T
        self.freq = freq
        self.lamb = 3e8/freq
        # Ground Track
        self.lat = lat
        self.lon = lon
        self.beam_radius = self.beam_radius_on_ground(self.Position.altitude, self.lamb, self.D_T)
        self.beam_area = self.beam_area_on_ground(self.beam_radius)
        self.beamwidth = self.effective_beam_width(self.beam_radius, self.Position.altitude)

    def beam_radius_on_ground(self, altitude, lamb, D_T):
        D_R = altitude * lamb/D_T
        beam_radius = D_R/2
        return beam_radius

    def beam_area_on_ground(self, rad_beam):
        return pi * rad_beam**2 

    def effective_beam_width(self, beam_radius, altitude):
        beamwidth = atan2(beam_radius, altitude)
        return beamwidth * 180/pi



def parseSamples(samples, lats, lons):
    i = 0
    for sample in samples:
        x = sample.x.value
        y = sample.y.value
        z = sample.z.value
        sat_pos = SphericalCoordinates.from_cartesian(x, y, z)
        yield SatPosEpoch(sat_pos, lats[i], lons[i])
        i = i + 1



def main(args):
    outputFile = args.output_file

    J2000_TT = time.Time("J2000", scale="tt")

    molniya = Orbit.from_classical(
        attractor=Earth,
        a=26570.5 * u.km,
        ecc=0.72 * u.one,
        inc=63.4 * u.deg,
        raan=0 * u.deg,
        argp=270 * u.deg,
        nu=0 * u.deg,
        epoch=J2000_TT
    )
    # The ground track provides us with 150 samples so we match that.
    num_samples = 360
    samples = molniya.sample(num_samples)
    period = molniya.period.value
    delta_t = period/num_samples
    s = ', '
    # Build spacecraft instance
    molniyaSpaceCraft = EarthSatellite(molniya, None)
    t_span = time_range(
        molniya.epoch, periods=num_samples, end=molniya.epoch + molniya.period.value * u.s
    )


    gp = GroundtrackPlotter()
    gp.update_layout(title="molniya ground track")

    # Plot previously defined EarthSatellite object
    gp.plot(
        molniyaSpaceCraft,
        t_span,
        label="Molniya GT",
        color="red",
        marker={
            "size": 10,
            "symbol": "triangle-right",
            "line": {"width": 1, "color": "black"},
        },
    )


    lats = gp.fig.data[1]['lat']
    lons = gp.fig.data[1]['lon']

    with open(outputFile, "w") as of:
        # for spacecraft velocity
        prev_r = None
        prev_theta = None
        prev_phi = None
        prev_lat = None
        prev_lon = None

        # for Ground Track 
        of.write(f"r{s}theta{s}phi{s}ground_track_beam_radius{s}ground_track_beam_area{s}v_r{s}v_theta{s}v_phi{s}beamwidth{s}v_groundtrack\n") #_theta{s}v_groundtrack_phi{s}v_relative_theta{s}v_relative_phi{s}v_lat{s}v_lon\n")
        i = 0
        for samp in parseSamples(samples, lats, lons):
            if prev_r is None or prev_theta is None or prev_phi is None:
                prev_r = samp.Position.r
                prev_theta = samp.Position.theta
                prev_phi = samp.Position.phi
                prev_lat = samp.lat
                prev_lon = samp.lon

            else:
                v_r = (samp.Position.r - prev_r)/delta_t
                v_theta = (samp.Position.theta - prev_theta)/delta_t
                v_phi = (samp.Position.phi - prev_phi)/delta_t
                # use lat and long to compute distance between cur samp and last samp
                v_groundtrack = geopy.distance.geodesic((samp.lat, samp.lon), (prev_lat, prev_lon)).m /delta_t
                of.write(f"{samp.Position.r}{s}"
                         f"{samp.Position.theta}{s}"
                         f"{samp.Position.phi}{s}"
                         f"{samp.beam_radius}{s}"
                         f"{samp.beam_area}{s}"
                         f"{v_r}{s}"
                         f"{v_theta}{s}"
                         f"{v_phi}{s}"
                         f"{samp.beamwidth}{s}"
                         f"{v_groundtrack}\n")
                prev_r = samp.Position.r
                prev_theta = samp.Position.theta
                prev_phi = samp.Position.phi
                prev_lat = samp.lat
                prev_lon = samp.lon

    
    # ORBIT
    myplt = molniya.plot()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                        prog = 'Space Ambient Power Satellite Properties Generator',
                        description = 'Builds csv file of data for orbital characteristics calcs',
                        epilog = 'Todo add some help menus')
    parser.add_argument("-f", "--frequency", help="transmit frequency ", default=5.8e9)
    parser.add_argument("-o", "--output_file", default='output.csv', help="output file name (csv format)")

    args = parser.parse_args()
    main(args)


