import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy import interpolate

# Define some units
micrometer = 1.0e-4 #cm
milimeter = 1.0e-1 #cm
centimeter = 1.0 #cm

# Class Material
class Material:
    def __init__(self, name="Material", density=1.0, stopping_power=lambda x: 1.0, file=False):
        self.name = name
        self.density = density #g/cm3
        if stopping_power == "file":
            stopping_power = self.load_stopping_power(file)
        self.stopping_power = np.vectorize(stopping_power) #MeV/(g/cm2)

    def load_stopping_power(self, file="./Material_Stopping_Power.py"):
        df = pd.read_csv(file, sep=';')
        self.df = df
        #print(df)
        energies = df["Kinetic_Energy"][1:].to_numpy()
        self.energies = energies.astype(np.float)
        total_stp_pow = df["Total.Stp.Pow."][1:].to_numpy()
        electronic_stp_pow = df["Electron.Stp.Pow."][1:].to_numpy()
        nuclear_stp_pow = df["Nuclear.Stp.Pow."][1:].to_numpy()
        
        self.total_stp_pow = total_stp_pow.astype(np.float) 
        self.electronic_stp_pow = electronic_stp_pow.astype(np.float) 
        self.nuclear_stp_pow = nuclear_stp_pow.astype(np.float)
        
        self.electronic_stopping_power = np.vectorize(interpolate.interp1d(self.energies, self.electronic_stp_pow, kind='linear')) 
        self.nuclear_stopping_power = np.vectorize(interpolate.interp1d(self.energies, self.nuclear_stp_pow, kind='linear')) 
        return interpolate.interp1d(self.energies, self.total_stp_pow, kind='linear')     



# Class Geometry
class GeometricalElement:
    def __init__(self, name="", seq_number=0, material=Material(), length=1.0, color='r'):
        self.name = name 
        self.seq_number = seq_number 
        self.material = material 
        self.length = length #cm
        self.dose = 0.0
        self.microdoses = [0.0]
        self.mean_dose = 0.0
        self.color = color
        self.npoints = 0
        self.initial_energy = 0.0
        self.final_energy = 0.0
        self.dx = 1.0e-4
        self.deposited_energies = [0.0]
        
    def restart():
        self.dose = 0.0

class Beam:
    def __init__(self, name="", energy=1.0, intensity=1.0, diameter=1.0):
        self.name = name
        self.energy = energy #MeV
        self.intensity = intensity #nA
        self.diameter = diameter #mm
        self.radius = self.diameter/2.0 #mm
        self.area = 1.0e-2*(np.pi*self.radius**2) #cm2
        
# Class Experiment
class Experiment:
    def __init__(self, name="", beam=Beam(), elements=[GeometricalElement()], irradiation_dose=5, irradiated_elements_names=[]):
        self.name = name
        self.beam = beam 
        self.elements = elements
        self.irradiation_dose = irradiation_dose
        self.irradiated_elements_names=irradiated_elements_names
        self.total_length = 0.0
    def add_element(self, element=GeometricalElement):
        self.elements.append(element)

    def compute(self, integration_density, time=1.0, verbose_level=0):
        print("Experiment Configuration:")
        print("Irradiation Dose:", self.irradiation_dose, "Gy")
        print("Beam Energy:", self.beam.energy, "MeV")
        print("Beam Current:", self.beam.intensity*1.0e3, "pA")
        print("Beam Diameter:", self.beam.diameter, "mm")
        print("Beam Area:", self.beam.area, "cm2")
        self.integration_density = integration_density
        initial_point = 0.0
        self.points = [0.0]
        for element in self.elements:
            element.initial_point = initial_point
            element.final_point = element.initial_point + element.length
            element.dx = 1.0/(element.material.density*self.integration_density)
            element.npoints = int(element.length/element.dx)
            print(element.name, element.length/micrometer, "micrometers", element.npoints, "points")
            self.total_length += element.length
            print(element.name, "dx = ", element.dx/micrometer, "micrometers")
            element.points = np.linspace(element.initial_point, element.final_point, element.npoints) 
            self.points = self.points + element.points.tolist()
            initial_point = element.final_point
        print("Total Length: ", self.total_length/micrometer, "micrometers")
        
        
        #self.points = np.linspace(0.0, self.dx*npoints, npoints+1)
        self.energy = self.beam.energy
        self.energies = [self.energy]
        self.deposited_energies = [0.0]
        self.deposited_energy_rates = [0.0]
        for element in self.elements:
            element.initial_energy = self.energy
            if self.energy <= 1.0e-3:
                self.energies.append(0.0)
                element.microdoses.append(0.0)
                break
            print(element.name, element.npoints, "points.")
            for i in range(element.npoints):
                let = element.material.stopping_power(self.energy)
                deposited_energy = let*element.material.density*(element.dx)
                self.deposited_energies.append(deposited_energy)
                element.deposited_energies.append(deposited_energy)
                self.deposited_energy_rates.append(deposited_energy/element.dx)
                self.energy = self.energy - deposited_energy
                if self.energy <= 1.0e-3:
                    element.final_energy = self.energy
                    self.energies.append(0.0)
                    element.microdoses.append((self.energy + deposited_energy)*(self.beam.intensity/self.beam.area)/(element.material.density*(element.dx)))
                    if verbose_level > 0.0:
                        print("Energy:", self.energy)
                        print("Deposited Energy:", deposited_energy)
                    break
                else:
                    self.energies.append(self.energy)
                    element.microdoses.append(deposited_energy*(self.beam.intensity/self.beam.area)/(element.material.density*(element.dx)))
                if verbose_level > 0.0:
                    print("Energy:", self.energy)
                    print("Deposited Energy:", deposited_energy)
            element.final_energy = self.energy
            element.dose = (self.beam.intensity/self.beam.area)*sum(element.deposited_energies)/(element.material.density*element.length)
            element.mean_dose = sum(element.microdoses)/element.npoints
           
              

    def show_results(self, shown_elements):
        while len(self.deposited_energies) < len(self.points):
            self.deposited_energies.append(0.0)
            self.energies.append(0.0)
        init_pos = 0.0 #self.elements[0].length
        total_points = 0
        for element in shown_elements:
            total_points = total_points + element.npoints
        for element in shown_elements:
            final_pos = init_pos + element.length
            plt.semilogy(self.points[-total_points:], self.deposited_energies[-total_points:], label="Deposited Energy")
            plt.axvspan(init_pos, final_pos, color=element.color, label=element.name, alpha=0.5)
            plt.semilogy(self.points[-total_points:], self.energies[-total_points:], label="Proton Energy")
            init_pos = final_pos
        print()
        print("---------------------------------------")
        for element in self.elements:
            print(element.name, "recibed a dose of", element.dose, "Gy per second")
            print(element.name, "recibed a mean dose of", element.mean_dose, "Gy per second")
            print(element.name, "first point recibed a microdose of", element.microdoses[1], "Gy per second")
            print(element.name, "last point recibed a microdose of", element.microdoses[-1], "Gy per second")
            print(element.name, "deposited energy per proton: ", sum(element.deposited_energies), "MeV")
        print()
        for element in self.elements:
            print(element.name, "Initial Proton Energy:", element.initial_energy, "MeV    ", "Final Proton Energy:", element.final_energy, "MeV")
            if element.name in self.irradiated_elements_names:
                if element.dose > 0.0:
                    self.irradiation_time = self.irradiation_dose/element.dose
                    self.irradiation_hours = int(self.irradiation_time/3600)
                    self.irradiation_minutes = int((self.irradiation_time-3600*self.irradiation_hours)/60)
                    self.irradiation_seconds = int(self.irradiation_time-3600*self.irradiation_hours-60*self.irradiation_minutes)
                    self.irradiation_miliseconds = int((self.irradiation_time-3600*self.irradiation_hours-60*self.irradiation_minutes-self.irradiation_seconds)*1.0e3)
                    print("Irradiation time required for a ", self.irradiation_dose, "Gy dose absorbed on", element.name, ":") 
                    print(self.irradiation_time, "s")
                    print(self.irradiation_hours, "h", self.irradiation_minutes, "min", self.irradiation_seconds, "s", self.irradiation_miliseconds, "ms")
                else:
                    print("No irradiation of ", element.name)
        print("---------------------------------------")
        plt.legend()
        plt.show()

    def save_results(self, outputfile):
        fig = plt.figure()
        while len(self.deposited_energies) < len(self.points):
            self.deposited_energies.append(0.0)
            self.energies.append(0.0)
        init_pos = 0.0
        for element in self.elements:
            plt.semilogy(self.points, self.deposited_energies, label="Deposited Energy")
            plt.axvspan(init_pos, init_pos + element.length, color='y', label=element.name + "Layer", alpha=0.5)
            plt.semilogy(self.points, self.energies, label="Proton Energy")
        print(file=open(outputfile, "a"))
        print("---------------------------------------", file=open(outputfile, "a"))
        print("Experiment:", self.name, file=open(outputfile, "a"))
        print(file=open(outputfile, "a"))
        print("Beam Energy:", self.beam.energy, file=open(outputfile, "a"))
        print("Beam Diameter:", self.beam.energy, file=open(outputfile, "a"))
        for element in self.elements:
            print(element.name, ":", element.length/micrometer, "micrometers", file=open(outputfile, "a"))
            print(element.name, "recibed a dose of", element.dose, "Gy per nA per second", file=open(outputfile, "a"))
        print("---------------------------------------", file=open(outputfile, "a"))
        plt.legend()
        fig_name = ".\DoseRates\ " + str(self.beam.energy) + "MeV.png"
        fig.savefig(fig_name)




