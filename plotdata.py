#import matplotlib.pyplot   as     plt
import csv
import argparse
import matplotlib.pyplot   as     plt


def main(args):

    t_vec = []
    r_vec = []
    v_gt_vec = []
    ratio_rad_gt_vec = []

    with open(args.filename, newline='') as fh:
        datareader = csv.reader(fh, delimiter=',')
        first = True;
        for row in datareader:
            if first:
                first = False
                continue
            try:
                r = float(row[0])
                t = float(row[-1])
                v_groundtrack = float(row[9])
                beam_radius = float(row[3])
            except ValueError:
                continue

            
            r_vec.append(r)
            t_vec.append(t)
            v_gt_vec.append(v_groundtrack)
            ratio_rad_gt_vec.append((beam_radius*2)/v_groundtrack)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    ax1.plot(t_vec, r_vec)
    ax1.set_ylabel("Radius of Orbit (m)")
    ax2.plot(t_vec, v_gt_vec)
    ax2.set_ylabel("Groundtrack Velocity (m/sec)")
    ax3.plot(t_vec, ratio_rad_gt_vec)
    ax3.set_ylabel("Beam Diameter to Groundtrack Velocity")
    ax3.set_xlabel("time (sec)")
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                        prog = 'Plot stuff',
                        description = 'plots csv stuff',
                        epilog = 'Todo add some help menus')
    parser.add_argument("-f", "--filename", help="csv file name")


    args = parser.parse_args()
    main(args)
