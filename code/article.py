# import necessary modules
import astropy.units as u
import huxt as H
import huxt_inputs as Hin
import numpy as np
import os
import h5py
import glob


def get_cr_number_list():
    """
    A function to compute the list of CRs that we have boundary conditon data for
    Returns:
    cr_num_list: a list of cr numbers (integers)
    """
    project_dirs = H._setup_dirs_()
    bc_path = os.path.join(project_dirs['boundary_conditions'], "*vr_r0.hdf")
    bc_files = glob.glob(bc_path)
    cr_num_list = []
    for f in bc_files:
        fn = os.path.basename(f)
        cr_tag = fn.split("_")[1]
        cr_num = int(cr_tag[2:])
        cr_num_list.append(cr_num)

    return cr_num_list


def build_cme_scenarios():
    """
    Build a dictionary of the average and fast CME scenarios as ConeCME objects.
    Returns:

    """
    # make CME scenarios
    cme_average = H.ConeCME(t_launch=0 * u.day, longitude=0.0 * u.deg, width=37.4 * u.deg, v=495 * (u.km / u.s),
                            thickness=0 * u.solRad)

    cme_fast = H.ConeCME(t_launch=0 * u.day, longitude=0.0 * u.deg, width=69.8 * u.deg, v=1070 * (u.km / u.s),
                         thickness=0 * u.solRad)

    cme_scenarios = {'cme_average': cme_average, 'cme_fast': cme_fast}
    return cme_scenarios


def run_experiment():
    # setup output file for data
    project_dirs = H._setup_dirs_()
    data_dir = project_dirs['output']
    out_path = os.path.join(data_dir, "CME_transit_data.hdf5")
    out_file = h5py.File(out_path, 'w')

    # specify carrington number range
    cr_num_list = get_cr_number_list()

    # Get ConeCME scenarios for experiment
    cme_scenarios = build_cme_scenarios()

    # Loop over carrington rotations
    for cr_number in cr_num_list:

        # Make group in output for CR
        cr_group = out_file.create_group('CR_{:d}'.format(cr_number))

        # import solar wind conditions
        vr_in = Hin.get_MAS_long_profile(cr_number, 0.0 * u.deg)

        # Export SW conditions and measure of variability in SW.
        cr_group.create_dataset('vr_in', data=vr_in)
        cr_group.create_dataset('vr_in_std', data=np.std(vr_in))

        # Loop over CME scenarios
        for cme_label, cme in cme_scenarios.items():

            cme_group = cr_group.create_group(cme_label)

            # create empty arrays to save variables
            lon_c = np.zeros(27)
            hit = np.zeros(27)
            hit_id = np.zeros(27)
            lon = np.zeros(27)
            r = np.zeros(27)
            t_arrive = np.zeros(27)
            t_transit = np.zeros(27)
            v = np.zeros(27)

            # Loop over daily longitudes in the carrington rotation
            for j in range(0, 27):

                # define dphi as the carrington longitude (value goes from 360 to 0 degrees)
                dphi = (360 * (27 - j) / 27) * u.deg

                # Initialise HUXt at this CR lon
                model = H.HUXt(v_boundary=vr_in, cr_num=cr_number, cr_lon_init=dphi, lon_start=-5 * u.deg,
                               lon_stop=5 * u.deg, simtime=7 * u.day, dt_scale=4)

                # Solve HUXt for this CME/CR lon
                model.solve([cme])
                cme_out = model.cmes[0]

                # Compute CME arrival at 1AU along nose lon of CME
                stats = cme_out.compute_arrival_at_location(0.0 * u.rad, 215.0 * u.solRad)

                # save outputs into arrays
                lon_c[j] = dphi.value
                hit[j] = stats['hit']
                hit_id[j] = stats['hit_id']
                lon[j] = stats['lon'].value
                r[j] = stats['r'].value
                t_arrive[j] = stats['t_arrive'].value
                t_transit[j] = stats['t_transit'].value
                v[j] = stats['v'].value

            # Save data to file
            out_data_dict = {'hit': hit, 'hit_id': hit_id, 'lon': lon, 'r': r, 't_arrive': t_arrive,
                             't_transit': t_transit, 'v': v, 'lon_c': lon_c}
            for key, val in out_data_dict.items():
                cme_group.create_dataset(key, data=val)

            out_file.flush()

    out_file.close()
    return


if __name__ == "__main__":
    run_experiment()
