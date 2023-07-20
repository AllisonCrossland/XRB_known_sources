import re
import pandas as pd
from datetime import datetime, timedelta

def dt_from_str(string):
    '''Extract datetime of swift zero epoch'''
    import re

    # Regular expression pattern to extract the date and time
    pattern = r'= (\d{4} \w{3} \d{2} at \d{2}:\d{2}:\d{2}.\d{3} UT)'

    # Search for the pattern in the string
    match = re.search(pattern, string)

    if match:
        # Extract the datetime substring
        datetime_string = match.group(1)

        # Define the format of the datetime substring
        format_string = '%Y %b %d at %H:%M:%S.%f UT'

        # Parse the datetime substring into a datetime object
        parsed_datetime = datetime.strptime(datetime_string, format_string)

#         # Print the parsed datetime
#         print(parsed_datetime)
#     else:
#         print("Datetime not found in the string.")
    return parsed_datetime


def dt_to_mjd(dt):
    # Calculate the Modified Julian Day (MJD)
    epoch = datetime(1858, 11, 17)  # MJD epoch
    mjd = (dt - epoch) / timedelta(days=1)

    # Return the Modified Julian Day (MJD)
    return mjd
    

def read_swift_data(fp):
    '''Read swift data file at filepath fp and returns photon counts and upper limits as dataframes'''
    
    # predefined values
    format_string = '%Y %b %d at %H:%M:%S.%f %Z'
    columns = ['Time', 'T_+ve', 'T_-ve', 'Rate', 'Ratepos', 'Rateneg', 'BGrate', 'BGerr', 'FracExp']
    
    # initialize variables
    WT_data = []
    PC_data = []
    upper_limits = []
    t0 = 0
    read_wt_data = 0
    read_pc_data = 0
    read_ul_data = 0
    
    # read field one line at a time
    with open(fp) as f:
        line = '1'
        while line:
            line = f.readline()
            
            # get zero point
            if 'Swift MET=' in line:
                t0 = dt_from_str(line)
            
            # end read at end of file (marked by Collapse)
            if 'Collapse' in line:
                line=False
                break
            
            # Append data to list if indicator variable is 1
            if read_wt_data:
                WT_data.append(line.split())
            
            if read_pc_data:
                PC_data.append(line.split())
                
            if read_ul_data:
                upper_limits.append(line.split())
                
            # Change window timing indicator to 1 after line with 'WT data'
            if '! WT data' in line:
                read_wt_data = 1
           
            # Change photon count indicator to 1 after line with 'PC data'
            if '! PC data' in line:
                read_pc_data = 1
                read_wt_data = 0

            # Change upper limit indicator to 1 after line with 'upper limits'
            if 'upper limits' in line:
                read_ul_data = 1
                read_wt_data = 0
                read_pc_data = 0
    
    # Convert saved lists to dataframes
    WT_data = pd.DataFrame(WT_data, columns=columns).dropna(subset=['FracExp'])
    PC_data = pd.DataFrame(PC_data, columns=columns).dropna(subset=['FracExp'])
    upper_limits = pd.DataFrame(upper_limits, columns=columns).dropna(subset=['FracExp'])
    
    # convert times to mjd
    WT_data['Time'] = dt_to_mjd(pd.to_timedelta(PC_data['Time'].astype('float'), unit='s') + t0)
    PC_data['Time'] = dt_to_mjd(pd.to_timedelta(PC_data['Time'].astype('float'), unit='s') + t0)
    upper_limits['Time'] = dt_to_mjd(pd.to_timedelta(upper_limits['Time'].astype('float'), unit='s') + t0)
    
    return PC_data, upper_limits, WT_data