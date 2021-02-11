import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy import interpolate

# Define some units
micrometer = 1.0e-4 #cm

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
        print(df)
        energies = df["Kinetic_Energy"][1:].to_numpy()
        self.energies = energies.astype(np.float)
        lets = df["Total.Stp.Pow."][1:].to_numpy()
        self.lets = lets.astype(np.float) 
        return interpolate.interp1d(self.energies, self.lets, kind='linear')     

# Define our materials:
paraffin = Material(name="Paraffin_Wax", density=0.9, stopping_power="file", file="./paraffin_proton_stopping_data.csv")
water = Material(name="Water", density=0.997, stopping_power="file", file="./water_proton_stopping_data.csv")


# Some plots
plt.loglog(water.energies, water.stopping_power(water.energies), label="Water LET")
plt.loglog(paraffin.energies, paraffin.stopping_power(paraffin.energies), "r:", label="Paraffin LET")
plt.legend()
plt.show()

my_energies = np.logspace(-1.0, 1.0, 100 )
plt.loglog(my_energies, water.stopping_power(my_energies), label="Water LET")
plt.loglog(my_energies, paraffin.stopping_power(my_energies), "r:", label="Paraffin LET")
plt.legend()
plt.show()

# Class Geometry
class GeometricalElement:
    def __init__(self, name="", seq_number=0, material=Material(), length=1.0):
        self.name = name 
        self.seq_number = seq_number 
        self.material = material 
        self.length = length #cm
        self.dose = 0.0
    
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
    def __init__(self, name="", beam=Beam(), elements=[GeometricalElement()]):
        self.name = name
        self.beam = beam 
        self.elements = elements
        self.total_length = 0.0
    def add_element(self, element=GeometricalElement):
        self.elements.append(element)

    def compute(self, npoints, time=1.0):
        for element in self.elements:
            print(element.name, element.length/micrometer, "micrometers")
            self.total_length += element.length
        print("Total Length: ", self.total_length/micrometer, "micrometers")
        self.dx = self.total_length / float(npoints) 
        print("dx = ", self.dx)
        self.points = np.linspace(0.0, self.dx*npoints, npoints)
        self.energy = self.beam.energy
        self.energies = [self.energy]
        self.deposited_energies = [0.0]
        for element in self.elements:
            if self.energy <= 1.0e-3:
                    break
            element_points = int(element.length/self.dx)
            print(element.name, element_points, "points.")
            for i in range(element_points):
                let = element.material.stopping_power(self.energy)
                deposited_energy = let*element.material.density*self.dx
                self.deposited_energies.append(deposited_energy)
                self.energy = self.energy - deposited_energy
                if self.energy <= 1.0e-3:
                    self.energies.append(0.0)
                    element.dose = element.dose + (self.energy + deposited_energy)*(self.beam.intensity/self.beam.area)
                    print("Energy:", self.energy)
                    print("Deposited Energy:", deposited_energy)
                    break
                else:
                    self.energies.append(self.energy)
                    element.dose = element.dose + deposited_energy*(self.beam.intensity/self.beam.area)
                print("Energy:", self.energy)
                print("Deposited Energy:", deposited_energy)
               

    def show_results(self):
        while len(self.deposited_energies) < len(self.points):
            self.deposited_energies.append(0.0)
            self.energies.append(0.0)
        plt.semilogy(self.points, self.deposited_energies, label="Deposited Energy")
        plt.axvspan(0.0, paraffin_layer.length, color='y', label="Paraffin Layer", alpha=0.5)
        plt.axvspan(paraffin_layer.length, paraffin_layer.length + cell_layer.length, color='m', label="Cell Layer", alpha=0.5)
        plt.semilogy(self.points, self.energies, label="Proton Energy")
        print()
        print("---------------------------------------")
        for element in self.elements:
            print(element.name, "recibed a dose of", element.dose, "Gy per nA per second")
        print("---------------------------------------")
        plt.legend()
        plt.show()

    def save_results(self, outputfile):
        fig = plt.figure()
        while len(self.deposited_energies) < len(self.points):
            self.deposited_energies.append(0.0)
            self.energies.append(0.0)
        plt.semilogy(self.points, self.deposited_energies, label="Deposited Energy")
        plt.axvspan(0.0, paraffin_layer.length, color='y', label="Paraffin Layer", alpha=0.5)
        plt.axvspan(paraffin_layer.length, paraffin_layer.length + cell_layer.length, color='m', label="Cell Layer", alpha=0.5)
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


# My Experiment
paraffin_layer = GeometricalElement(name="Paraffin", length=20.0*micrometer, material=paraffin)
cell_layer = GeometricalElement(name="Cells", length=10.0*micrometer, material=water)
my_beam = Beam(energy=10.0, intensity=1.0, diameter=6.0)
my_experiment = Experiment(name="Paraffin and Cells", beam=my_beam, elements=[paraffin_layer, cell_layer])
my_experiment.compute(10000)
my_experiment.show_results()
my_experiment.save_results(outputfile="./results.txt")

