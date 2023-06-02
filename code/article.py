# import necessary modules
import astropy.units as u
import huxt as H
import huxt_inputs as Hin
import csv
import numpy as np
import os

def run_experiment():
    
    # specify carrington number range
    cr_min = 1625
    cr_max = 2265

    # create a loop to simulate solar wind conditions in each carrington number
    for i in range(cr_min, cr_max):
    
    
        # create empty arrays to save variables
        num_c=np.zeros(27)
        lon_c=np.zeros(27)

        hit=np.zeros(27)
        hit_id=np.zeros(27)
        lon=np.zeros(27)
        r=np.zeros(27)
        t_arrive=np.zeros(27)
        t_transit=np.zeros(27)
        v=np.zeros(27)
        
        # calculate transit time if solar wind conditions are available
        try:
            
            # import solar wind conditions
            vr_in = Hin.get_MAS_long_profile(i, 0.0*u.deg)
            
            # create a loop for each longitude in the carrington rotation
            for j in range(0, 27):
                
                # define dphi as the carrington longitude (value goes from 360 to 0 degrees)
                dphi = (360*(27-j)/27)*u.deg
                
                # simulate ambient solar wind conditions
                model = H.HUXt(v_boundary=vr_in, cr_num=i, cr_lon_init=dphi, lon_start=-5*u.deg,lon_stop=5*u.deg,simtime=7*u.day, dt_scale=4)
            
                # introduce a spherical CME
                cme = H.ConeCME(t_launch=0*u.day, longitude=0.0*u.deg, width=37.4*u.deg, v=495*(u.km/u.s), thickness=0*u.solRad)
                cme_list = [cme]
                model.solve(cme_list)
                cme = model.cmes[0]
                stats = cme.compute_arrival_at_location(0.0*u.rad, 215.0*u.solRad)
                
                # save outputs into arrays
                num_c[j]=i
                lon_c[j]=dphi.value
                
                hit[j]=stats['hit']
                hit_id[j]=stats['hit_id']
                lon[j]=stats['lon'].value
                r[j]=stats['r'].value
                t_arrive[j]=stats['t_arrive'].value
                t_transit[j]=stats['t_transit'].value
                v[j]=stats['v'].value
            
        # otherwise save transit time values as nan value
        except:
            
            # saving outputs with nan values in transit time                
            num_c[:]=i
            hit = hit*np.nan
            hit_id = hit_id*np.nan
            lon = lon*np.nan
            r = r*np.nan
            t_arrive = t_arriva*np.nan
            t_transit = t_transit*np.nan
            v = v*np.nan
            lon_c = 360*(27 - np.arange(0,27))/27
                
                
        # save the data into a csv        
        finally:
            
            # Save data to file
            project_dirs = H._setup_dirs_()
            data_dir = project_dirs['output']
            out_path = os.path.join(data_dir,"cr_%s.csv")
            # open a new csv file in a specified path
            with open(out_path%i, 'w', newline='') as df:
                writer=csv.writer(df)
                
                # write in the column names
                writer.writerow(('num_c','lon_c','hit','hit_id','lon','r','t_arrive','t_transit','v'))
                    
                # write in the values
                for k in range(0,27):
                    writer.writerow((num_c[k],lon_c[k],hit[k],hit_id[k],lon[k],r[k],t_arrive[k],t_transit[k],v[k])) 
                    
        break
    return
    
if __name__ == "__main__":
    run_experiment()