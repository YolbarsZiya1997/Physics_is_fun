BOHR_TO_A = 0.5291772
HARTREE_TO_EV = 27.211396132
 
class Cube:
    def __init__(self, filename, fermi_level, band_gap):
        self.filename = filename
        self.name = filename.split('.')[0]
        self.valence_band = fermi_level
        self.conduction_band = fermi_level + band_gap
         
        if filename.split('.')[-1] == 'cube':
            self.parse_cube_file()
            self.project_along_z()
            self.write_projection_file()
         
        elif filename.split('.')[-1] == 'dat':
            self.parse_projection_file()
         
        else:
            error_message = 'File format of file {:s} not supported'.format(self.filename)
            raise Exception(error_message)
             
    def parse_cube_file(self):
        with open(self.filename, 'r') as f:
            raw_data = f.readlines()
             
        self.natoms = int(raw_data[2].split()[0])
        self.nx, self.ny, self.nz = [int(line.split()[0]) for line in raw_data[3:6]]
        self.dx, self.dy, self.dz = [float(line.split()[i + 1]) for i, line in enumerate(raw_data[3:6])]
        self.A = self.nx * self.dx * BOHR_TO_A
        self.B = self.ny * self.dy * BOHR_TO_A
        self.C = self.nz * self.dz * BOHR_TO_A
         
        self.data = []
        for line in raw_data[self.natoms + 6:]:
            self.data.extend(float(x) for x in line.split())
             
        self.potential = sp.zeros((self.nx, self.ny, self.nz))
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    self.potential[i, j, k] = self.data[i * self.ny * self.nz + j * self.nz + k] * HARTREE_TO_EV
         
        print('Cube file {:s} parsed\n'.format(self.filename))
 
    def parse_projection_file(self):
        with open(self.filename, 'r') as f:
            raw_data = f.readlines()
             
        self.natoms = int(raw_data[0].split()[0])
        self.nx, self.ny, self.nz = [int(value) for value in raw_data[1].split()]
        self.dx, self.dy, self.dz = [float(value) for value in raw_data[2].split()]
        self.A = self.nx * self.dx * BOHR_TO_A
        self.B = self.ny * self.dy * BOHR_TO_A
        self.C = self.nz * self.dz * BOHR_TO_A
         
        self.projection_z = sp.zeros(self.nz)
        for i, line in enumerate(raw_data[3:]):
            self.projection_z[i] = float(line.split()[0])
         
        print('Cube file {:s} parsed\n'.format(self.filename))
         
    def write_projection_file(self):
        with open(self.name + '.dat', 'w') as f:
            f.write('{:d}\n'.format(self.natoms))
            f.write('{:d} {:d} {:d}\n'.format(self.nx, self.ny, self.nz))
            f.write('{: 10.7f} {: 10.7f} {: 10.7f}\n'.format(self.dx, self.dy, self.dz))
            for value in self.projection_z:
                f.write('{: 8.5f}\n'.format(value))
         
        print('Projection file {:s} written\n'.format(self.name + '.dat'))
     
    def project_along_z(self):
        self.projection_z = sp.sum(sp.sum(self.potential, axis=0), axis=0)
        self.projection_z /= (self.nx * self.ny)
     
    def plot_projection(self):
        x_axis = sp.linspace(0, self.C, self.nz)
        plt.plot(x_axis, self.projection_z, label=self.name)
        plt.title(self.name)
        plt.xlim([min(x_axis), max(x_axis)])
        plt.ylim([min(self.projection_z) - 5, max(self.projection_z) + 5])
        plt.show()
         
class Cube_Interface(Cube):
    def __init__(self, filename, fermi_level, band_gap, start_peak=1, end_peak=3, cutoff=-10.0):
        super().__init__(filename, fermi_level, band_gap)
        self.start_peak = start_peak
        self.end_peak = end_peak
        self.cutoff = cutoff
     
    def find_peaks(self):
        peaks = []
        last_value = self.projection_z[0]
        is_looking_for_minimum = True
         
        for i, value in enumerate(self.projection_z):
            if not is_looking_for_minimum and value > self.cutoff:
                is_looking_for_minimum = True
                 
            if is_looking_for_minimum:
                if value > last_value and value < self.cutoff:
                    peaks.append(i - 1)
                    is_looking_for_minimum = False
                     
                last_value = value
             
        return peaks
     
    def calculate_average_potential(self):
        peaks = self.find_peaks()
        center_potential = self.projection_z[peaks[self.start_peak] : peaks[self.end_peak]]
        self.average_potential = sum(center_potential) / len(center_potential)
        self.vacuum_level = max(self.projection_z)
         
    def plot_projection(self):
        peaks = self.find_peaks()
        start_peak_pos = peaks[self.start_peak] * self.dz * BOHR_TO_A
        end_peak_pos = peaks[self.end_peak] * self.dz * BOHR_TO_A
        plt.hlines(self.average_potential, start_peak_pos, end_peak_pos)
        super().plot_projection()
         
class Cube_Bulk(Cube):
    def __init__(self, filename, fermi_level, band_gap):
        super().__init__(filename, fermi_level, band_gap)
     
    def calculate_average_potential(self):
        self.average_potential = sum(self.projection_z) / len(self.projection_z)
     
    def plot_projection(self):
        plt.hlines(self.average_potential, 0.0, self.C)
        plt.hlines(self.valence_band, 0.0, self.C)
        plt.hlines(self.conduction_band, 0.0, self.C)
        super().plot_projection()
         
class Alignment:
    def __init__(self, bulk, interface):
        self.bulk = bulk
        self.interface = interface
         
        self.name = bulk.name.split('_')[0]
         
        if not self.interface.average_potential:
            self.interface.calculate_average_potential()
             
        if not self.bulk.average_potential:
            self.bulk.calculate_average_potential()
         
        self.calculate_band_edges()
         
    def calculate_band_edges(self):
        vacuum = self.interface.vacuum_level
        valence_band = self.bulk.valence_band - self.bulk.average_potential
        conduction_band = self.bulk.conduction_band - self.bulk.average_potential
        self.adjusted_valence_band = (self.interface.average_potential + valence_band) - vacuum
        self.adjusted_conduction_band = (self.interface.average_potential + conduction_band) - vacuum
        print('Conduction band: {: 5.3f}'.format(self.adjusted_conduction_band))
        print('Valence band:    {: 5.3f}'.format(self.adjusted_valence_band))
     
 
TiO2_vacuum = Cube_Interface('TiO2_vacuum.dat', -2.655171, 4.606532)
TiO2_vacuum.calculate_average_potential()
TiO2_vacuum.plot_projection()
TiO2_bulk = Cube_Bulk('TiO2_bulk.dat', 4.601682, 4.674498)
TiO2_bulk.calculate_average_potential()
TiO2_bulk.plot_projection()
TiO2 = Alignment(TiO2_bulk, TiO2_vacuum)