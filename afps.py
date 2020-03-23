import csv
import os
from collections import namedtuple
import matplotlib.pyplot as plt

import time
def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print ('%r  %2.2f ms' % \
                  (method.__name__, (te - ts) * 1000))
        return result

    return timed

class Synthesis:
    "Use the .pep data folder to import your data into the Synthesis object."
    def __init__(self,pepfile):
        self.pepfile = pepfile
        self.serial_no, self.sequence, self.datafile  = self.get_afps_sequence(pepfile)
        self.afps_data = self.parse_afps_datafile()
        self.integrals = self.get_integrals()

    def print_integrals(self):
        for number,amino,area,height,width in self.integrals:
            print(number, amino, area, height, width)

    @timeit
    def parse_afps_datafile(self):
        "Parse the afps datafile, returning a list of dictionaries containing each timepont's data"
        with open(self.datafile) as csvfile:
            csvreader = csv.reader(csvfile)
            
            next(csvreader,None)
            sequence = ''.join(next(csvreader)).split(": ")[1]

            for i in range(3):
                print(next(csvreader,None))
            
            headers = next(csvreader)
            # Put the data in a dictionary. date and time values should be treated as strings,
            # step_no as an int, and the rest as floats. 
            afps_data = [ {key:value if key in ['DATE_YYYYMMD', 'TIME_HHMMSSS.SSS'] \
                                     else int(value.split('.')[0]) if key == 'STEP_NO'\
                                     else float(value) for (key,value) in zip(headers,row)} \
                                     for row in csvreader ]
	
            for row in afps_data:
                hrs = int(row['TIME_HHMMSSS.SSS'][0:2])
                mins = int(row['TIME_HHMMSSS.SSS'][2:4])
                secs = int(row['TIME_HHMMSSS.SSS'][4:6])
                ms = int(row['TIME_HHMMSSS.SSS'][7:10])
                
                float_ms = float(hrs*60*60*1000 + mins*60*1000 + secs*1000 + ms)
                
                row['TIME_MS'] = float_ms
                print(float_ms)
            
            return afps_data
    
    @staticmethod
    def get_afps_sequence(pepfile):
        "Returns (sequence, serial_no, datafile) for a .pep file input."

        if not pepfile.endswith('/'): pepfile += '/'
        AFPS_File = namedtuple('AFPS_File', 'serial_no sequence afps_datafile_path')
        for afps_file in os.listdir(pepfile):
            if afps_file.split('.')[-1] == 'afps':
                # TODO: Make SN an integer
                serial_no = int(afps_file.split('.')[0].split('_')[-1].strip('SN0'))
                afps_datafile_path = pepfile + afps_file
                
                with open(afps_datafile_path) as csvfile:
                    csvreader = csv.reader(csvfile)
                    next(csvreader,None)
                    sequence = ''.join(next(csvreader)).split(": ")[1]
                
                return AFPS_File(serial_no = serial_no, \
                                 sequence = sequence,   \
                                 afps_datafile_path = afps_datafile_path)

    @staticmethod
    def get_pep_files(path):
        "Yields .pep files in path." 
        if not path.endswith('/'): path = path + '/'
        for afps_dir in os.listdir(path):
            if afps_dir.split('.')[-1] == 'pep':
                yield(path + afps_dir)

    @classmethod
    def from_serial(cls, serial_no, path='./'):
        "Create a Synthesis object from serial number."
        if not path.endswith('/'): path += '/'
        for item in os.listdir(path):
            extension = item.split('.')[-1]
            name = item.split('.')[0]
            if extension == 'pep':
                if int(name.split('_')[-1].strip('SN')) == serial_no:
                    print('hello')
                    return cls(path+item)
        print('Serial Number', serial_no, 'not found')
        return None

    def plot_uv(self):
        "Plots the raw UV data using matplotlib"
        
        time,uv,step = zip(*((row['TIME_MS'], row['UV_VIS'], row['STEP_NO']) for row in self.afps_data))
        fig, line = plt.subplots()
        line.plot(time,uv)
        plt.show()
        
    def plot_pressure(self):
        "Plots the raw pressure data using matplotlib"
        
        time,pressure,step = zip(*((row['TIME_MS'], row['PRESSURE_1'], row['STEP_NO']) for row in self.afps_data))
        fig, line = plt.subplots()
        line.plot(time,pressure)
        plt.show()

    def plot_temp(self):
        "Plots the raw temperature data using matplotlib"
        
        time,temp4,temp5,temp6,temp8,step = zip(*((row['TIME_MS'],row['TC_4'],row['TC_5'],row['TC_6'], row['TC_8'], row['STEP_NO']) for row in self.afps_data))
        fig, line = plt.subplots()
        line.plot(time,temp4, label="Zone 4")
        line.plot(time,temp5, label="Zone 5")
        line.plot(time,temp6, label="Zone 6")
        line.plot(time,temp8, label="Zone 8")
        line.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
        plt.show()

    def plot(self,field):
        time,field,step =  zip(*((row['TIME_MS'], row[field], row['STEP_NO']) for row in self.afps_data))
        fig, line = plt.subplots()
        line.plot(time,field)
        plt.show()
    
    def plot_step(self,field,step):
        time,field,step =  zip(*((row['TIME_MS'], row[field], row['STEP_NO']) for row in self.afps_data if row['STEP_NO'] == step))
        fig, line = plt.subplots()
        line.plot(time,field)
        plt.show()


    def get_data(self,*fields):
        # Returns the specified columns
        return tuple(zip(*([row[field] for field in fields] for row in self.afps_data)))

    def correct_uv_baseline(self, plot=False):
        "Returns a baseline corrected UV trace. Requires peakutils and numpy."
        import peakutils
        import numpy as np

        time, step_no, numpy_uv = np.array(self.get_data('TIME_MS', 'STEP_NO', 'UV_VIS'))
        
        baseline = peakutils.baseline(numpy_uv,3)
        numpy_uv_corrected = numpy_uv - baseline

        return list(zip(step_no,time,numpy_uv_corrected))

    def detect_depro_steps(self):
        # Parse out UV peaks by step number. 42, 64, 88....
        all_depro_data = [ row for row in self.afps_data           \
                               if ( row['STEP_NO'] - 14 )% 22 == 0 \
                               and row['STEP_NO'] >= 36 ]

        # Create a set of all unique step numbers captured here (one for each deprotection peak)
        depro_steps = {row['STEP_NO'] for row in all_depro_data}


        # Sorry for all of the comprehensions. 
        # This puts all of the data belonging to each depro peak
        # in a dictionary indexed by the step number it belongs to.
        #depro_peaks = {step_no: [rows for rows in all_depro_data        \
        #                              if rows['STEP_NO'] == step_no]    \
        #                              for step_no in depro_steps}

        # How many steps did we get?
        print(len(depro_steps), 'deprotection peaks found')# How many residues in our sequence?
        print(len(self.sequence), 'residues in sequence')

        return depro_steps

    def integrate(self, peak):
        # Performs the riemann sum underneath the given curve with format [*[time,data]]. Returns area, FWHM, max.
        # Sort the data by time
        peak.sort(key=lambda x: x[0])
        
        # Iterate over the list and compute the sum
        sum = 0
        for index,item in enumerate(peak):
            if index == 0:
                pass
            else:
                dt = peak[index][0] - peak[index - 1][0]
                avg_uv = (peak[index][1] + peak[index - 1][1])/2
                sum += dt * avg_uv
        return sum
     
    def get_max(self,peak):
        time,uv = zip(*peak)
        
        maximum = max(uv)
        half_max = maximum / 2

        for time,uv in peak:
            if uv - half_max > 0:
                left = time

        for time,uv in reversed(peak):
            if uv - half_max > 0:
                right = time
        # time is in ms; return in units of seconds
        width = abs(right - left) / 1000
        return maximum, width

    def get_integrals(self):
        steps = self.detect_depro_steps()
        corrected_uv_data = self.correct_uv_baseline()
        integrals = []

        for index,step in enumerate(sorted(steps)):
            peak = [ [i[1], i[2]] for i in corrected_uv_data if i[0] == step ]
            integral = self.integrate(peak)
            maximum,width = self.get_max(peak)
            integrals.append((index,self.sequence[-1 - index],integral,maximum,width))
        
        return integrals
