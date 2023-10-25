# import necessary modules
import astropy.units as u
import huxt as H
import huxt_inputs as Hin
import numpy as np
import os
import h5py


def run_experiment():
    # setup output file for data
    project_dirs = H._setup_dirs_()
    data_dir = project_dirs['output']
    out_path = os.path.join(data_dir, "CME_transit_data.hdf5")
    out_file = h5py.File(out_path, 'w')

    # specify carrington number range
    cr_min = 1625
    cr_max = 2265

    # create a loop to simulate solar wind conditions in each carrington number
    for cr_number in range(cr_min, cr_max):

        # create empty arrays to save variables
        num_c = np.zeros(27)
        lon_c = np.zeros(27)

        hit = np.zeros(27)
        hit_id = np.zeros(27)
        lon = np.zeros(27)
        r = np.zeros(27)
        t_arrive = np.zeros(27)
        t_transit = np.zeros(27)
        v = np.zeros(27)

        # calculate transit time if solar wind conditions are available
        try:

            # import solar wind conditions
            vr_in = Hin.get_MAS_long_profile(cr_number, 0.0 * u.deg)

            # create a loop for each longitude in the carrington rotation
            for j in range(0, 27):
                # define dphi as the carrington longitude (value goes from 360 to 0 degrees)
                dphi = (360 * (27 - j) / 27) * u.deg

                # simulate ambient solar wind conditions
                model = H.HUXt(v_boundary=vr_in, cr_num=cr_number, cr_lon_init=dphi, lon_start=-5 * u.deg,
                               lon_stop=5 * u.deg, simtime=7 * u.day, dt_scale=4)

                # introduce a climatological average CME
                cme = H.ConeCME(t_launch=0 * u.day, longitude=0.0 * u.deg, width=37.4 * u.deg, v=495 * (u.km / u.s),
                                thickness=0 * u.solRad)
                cme_list = [cme]
                model.solve(cme_list)
                cme = model.cmes[0]
                stats = cme.compute_arrival_at_location(0.0 * u.rad, 215.0 * u.solRad)

                # save outputs into arrays
                num_c[j] = cr_number
                lon_c[j] = dphi.value

                hit[j] = stats['hit']
                hit_id[j] = stats['hit_id']
                lon[j] = stats['lon'].value
                r[j] = stats['r'].value
                t_arrive[j] = stats['t_arrive'].value
                t_transit[j] = stats['t_transit'].value
                v[j] = stats['v'].value

        # otherwise save transit time values as nan value
        except:

            # saving outputs with nan values in transit time                
            num_c[:] = cr_number
            hit = hit * np.nan
            hit_id = hit_id * np.nan
            lon = lon * np.nan
            r = r * np.nan
            t_arrive = t_arrive * np.nan
            t_transit = t_transit * np.nan
            v = v * np.nan
            lon_c = 360 * (27 - np.arange(0, 27)) / 27

        # save the data into a csv        
        finally:

            # Save data to file
            cr_group = out_file.create_group('CR_{:d}'.format(cr_number))
            out_data_dict = {'num_c': num_c, 'hit': hit, 'hit_id': hit_id, 'lon': lon, 'r': r, 't_arrive': t_arrive,
                             't_transit': t_transit, 'v': v, 'lon_c': lon_c}

            for key, val in out_data_dict.items():
                cr_group.create_dataset(key, data=val)

            out_file.flush()

    out_file.close()
    return


if __name__ == "__main__":
    run_experiment()
