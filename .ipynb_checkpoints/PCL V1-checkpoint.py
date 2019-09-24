from PyQt4.QtCore import *
from PyQt4.QtGui import *
import sys
import numpy as np
import random
import networkx as nx
import community
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector
import numpy as np
from skimage.data import coins
from os import listdir
from scipy import stats
from skimage import io


sys.path.append(r"C:\Users\matth\OneDrive\Documents\Python Code\Matts Modules")
import matts_admin_functions
import matts_network_functions
import matts_matrix_functions
import PCL_Parameters

#Global Variables
global number_of_timepoints
global number_of_neurons

#Thresholds
global cell_probability_threshold
global cell_noise_threshold



#Neuron Data
global region_images
global tiff_dimensions
global neuron_coordinates
global flourescence_traces
global cell_probabilities
global cell_noises
global spike_matricies
global calcium_traces
global baseline_traces

#Decision variables
global cell_probability_decisions
global within_region_decisions
global cell_noise_decisions
global neuron_decisions
global region_assignments

#Lasso Variables
global lassos
global region_combo_boxes
global currently_selected_region

#GUI Variables
global cell_trace_timepoint
global current_cell_trace_cell
global current_cell_trace_plane
global current_cell_trace_file
global current_cell_outline_x_coords
global current_cell_outline_y_coords
global window_radius

global cell_trace_figures
global cell_trace_canvases
global position_indicator_axes


#BCL Parameters
global wdt
global lamb
global varB
global varC
global initial_baseline
global Cmean
global fourier_cutoff


"""
Tectum Parameters
Spiking probability 0.2
Calcium Decay 0.5
Baseline Variance 0.01
Calcium Variance 1
Initial baseline 
Mean calcium 
Fourier Cutoff 20000
frequency 4.8
"""


"""
Whole Brain Parameters
Spiking probability 0.2
Calcium Decay 2
Baseline Variance 0.01
Calcium Variance 0.5
Initial baseline 1
Mean calcium 1
Fourier Cutoff 1000
frequency 2.7

"""

#BCL Parameters
wdt = 0.1
lamb = 0.5
varB = 0.013#0.01
varC = 1
initial_baseline = 0
Cmean = 0
frequency = 2.7
fourier_cutoff = 8500


window_radius = 10
image_radius = 25
cell_trace_timepoint = 0

folder_location = r"F:\F1_tct_suite2p"
folder_location = r"C:\Users\matth\OneDrive\Documents\New Tectal Data\5p_tectum_suite2p"
output_file_location = r"C:\Users\matth\OneDrive\Documents\output.csv"
coordinate_output_file_location = r"C:\Users\matth\OneDrive\Documents\active_rois.csv"

number_of_planes = 1


cell_probability_threshold = 0.8


def write_list_to_file(list, filename):
    with open(filename, "w") as outfile:
        for entries in list:
            outfile.write(entries)
            outfile.write("\n")

def calculate_number_of_cells():
    number_of_cells = 0
    for plane in range(number_of_planes):
        for cell in range(len(neuron_decisions[plane])):
            if neuron_decisions[plane][cell] == True:
                number_of_cells += 1

    main.number_of_neurons_label.setText("Number Of Neurons: " + str(number_of_cells))


def threshold_cells():
    global neuron_decisions
    for plane in range(number_of_planes):
        for cell in range(len(neuron_decisions[plane])):
            if cell_probability_decisions[plane][cell] and within_region_decisions[plane][cell] and cell_noise_decisions[plane][cell]:
                neuron_decisions[plane][cell] = True
            else:
                neuron_decisions[plane][cell] = False



def setup_variables():
    global neuron_coordinates
    global neuron_decisions
    global within_region_decisions
    global cell_probability_decisions
    global cell_noises
    global cell_noise_decisions
    global region_assignments
    global full_neuron_regions
    global spike_matricies
    global cell_trace_figures
    global position_indicator_axes
    global cell_trace_canvases
    global calcium_traces
    global baseline_traces
    global region_combo_boxes
    global currently_selected_region
    global neuron_locations
    global fourier_traces

    neuron_coordinates          = matts_admin_functions.create_empty_list(number_of_planes)
    neuron_decisions            = matts_admin_functions.create_empty_list(number_of_planes)
    within_region_decisions     = matts_admin_functions.create_empty_list(number_of_planes)
    cell_probability_decisions  = matts_admin_functions.create_empty_list(number_of_planes)
    cell_noises                 = matts_admin_functions.create_empty_list(number_of_planes)
    cell_noise_decisions        = matts_admin_functions.create_empty_list(number_of_planes)
    region_assignments          = matts_admin_functions.create_empty_list(number_of_planes)
    spike_matricies             = matts_admin_functions.create_empty_list(number_of_planes)
    cell_trace_figures          = matts_admin_functions.create_empty_list(number_of_planes)
    cell_trace_canvases         = matts_admin_functions.create_empty_list(number_of_planes)
    position_indicator_axes     = matts_admin_functions.create_empty_list(number_of_planes)
    calcium_traces              = matts_admin_functions.create_empty_list(number_of_planes)
    baseline_traces             = matts_admin_functions.create_empty_list(number_of_planes)
    region_combo_boxes          = matts_admin_functions.create_empty_list(number_of_planes)
    neuron_locations            = matts_admin_functions.create_empty_list(number_of_planes)
    fourier_traces              = matts_admin_functions.create_empty_list(number_of_planes)
    currently_selected_region   = matts_admin_functions.create_list_of_item(number_of_planes,"Forebrain")


    print "Number of neurons", number_of_neurons

    for plane in range(number_of_planes):
        for neuron in range(number_of_neurons[plane]):
            neuron_decisions            [plane].append(True)
            within_region_decisions     [plane].append(False) #change back later
            cell_probability_decisions  [plane].append(True)
            cell_noises                 [plane].append(0)
            cell_noise_decisions        [plane].append(True)
            region_assignments          [plane].append("N/A")
            spike_matricies             [plane].append([])
            cell_trace_figures          [plane].append([])
            position_indicator_axes     [plane].append([])
            cell_trace_canvases         [plane].append([])
            calcium_traces              [plane].append([])
            baseline_traces             [plane].append([])
            fourier_traces              [plane].append([])
            neuron_locations            [plane].append(0)

    populate_trace_combobox()



def load_data():
    global lassos
    global folder_location
    global number_of_neurons

    #Get Location Of Suite2P Data
    folder_location = QFileDialog.getExistingDirectory(main, "Open a folder")

    #Determine The Number Of Neurons
    get_neuron_numbers()

    #Setup Variables Ready To Be Populated
    setup_variables()

    #Populate GUI
    populate_gui()

    #Load The Mean Images
    load_mean_images()

    #Load The Neuron Coordinates
    load_rois()

    #Load The Cell Probabilities
    load_cell_probabilities()

    #Load The Raw Flouresence Traces
    load_flourescence_traces()


    #Average Flourescence Traces
    #average_flourescence_traces()

    #Create The Region Selection Lassos
    lassos = create_lassos()

    #Draw All The Neurons
    draw_all_roi_scenes()


def average_trace(trace,window):
    average_trace = []

    for timepoint in range(len(trace)):

        if timepoint + window == len(trace):
            upper_bound = (len(trace) - timepoint)
        else:
            upper_bound = window

        average = np.mean(trace[timepoint : timepoint+upper_bound])
        average_trace.append(average)

    return average_trace




def get_neuron_numbers():
    global number_of_neurons
    number_of_neurons = []

    for plane in range(number_of_planes):
        stats_file_location = str(folder_location + str(r"\plane"+str(plane)) + "\stat.npy")
        stats = np.load(stats_file_location)
        plane_size = len(stats)
        number_of_neurons.append(plane_size)


def populate_gui():
    main.populate_region_selection_tab()



def load_rois():
    global neuron_coordinates

    for plane in range(number_of_planes):
        stats_file_location = str(folder_location + str(r"\plane"+str(plane)) + "\stat.npy")
        stats = np.load(stats_file_location)

        for neuron in range(number_of_neurons[plane]):
            neuron_coordinates[plane].append([stats[neuron]["med"][0],stats[neuron]["med"][1],plane])









def populate_trace_combobox():

    for plane in range(number_of_planes):
        main.traces_combox_box.addItem("Plane " + str(plane))

    max = len(neuron_decisions[0])
    main.neuron_end_spinner.setMaximum(max)
    main.neuron_start_spinner.setMaximum(max)

def load_cell_probabilities():
    global cell_probabilities
    cell_probabilities = matts_admin_functions.create_empty_list(number_of_planes)

    for plane in range(number_of_planes):
        is_cell_file_location = str(folder_location + str(r"\plane"+str(plane)) +"\iscell.npy")
        is_cell = np.load(is_cell_file_location)

        for cell in range(len(is_cell)):
            cell_probabilities[plane].append(is_cell[cell][1])


def calculate_highest_minimums():
    global cell_noises
    sliding_window_size = 9

    for plane in range(number_of_planes):

        main.minimum_values_progress_bar.setValue((float(plane+1) / number_of_planes) * 100)
        app.processEvents()

        for cell in range(len(flourescence_traces[plane])):
            if neuron_decisions[plane][cell] == True:
                stop_point = len(flourescence_traces[plane][cell])-sliding_window_size
                minimum_values = []

                for timepoint in range(0, stop_point,sliding_window_size):
                    timeseries = flourescence_traces[plane][cell][timepoint:timepoint + sliding_window_size]
                    minimum = np.min(timeseries)
                    minimum_values.append(minimum)

                highest_minimum = np.max(minimum_values)
                cell_noises[plane][cell] = highest_minimum


    plot_highest_minimums()

def plot_highest_minimums():
    highest_minimum_list = []
    for plane in range(number_of_planes):
        for cell in range(len(cell_noises[plane])):
            if neuron_decisions[plane][cell] == True:
                highest_minimum_list.append(cell_noises[plane][cell])

    plot_histogram(main.minimum_values_figure, main.minimum_values_canvas, highest_minimum_list, "Highest Minimum Values")





def load_flourescence_traces():
    global flourescence_traces
    global number_of_timepoints
    flourescence_traces = matts_admin_functions.create_empty_list(number_of_planes)

    for plane in range(number_of_planes):
        trace_file_location = str(folder_location + str(r"\plane"+str(plane)) + "\F.npy")
        traces = np.load(trace_file_location)
        flourescence_traces[plane] = traces

    number_of_timepoints = len(flourescence_traces[plane][0])




def update_cell_trace_tiff_index(image,coordinates):
    global current_cell_trace_cell
    global current_cell_trace_plane
    global cell_trace_figures
    global cell_trace_canvases
    global position_indicator_axes


    tiff_index = cell_trace_timepoint % 200

    window_open = tiff_index - window_radius
    if window_open < 0:
        window_open = 0

    window_close = tiff_index + window_radius
    if window_close >= 200:
        window_close = 199

    image_matrix = image[window_open:window_close]
    image_matrix = np.array(image_matrix)
    image_matrix = np.mean(image_matrix, axis=0)
    image_matrix = matts_matrix_functions.extract_matrix_region_around_point(image_matrix,coordinates,image_radius)
    image_matrix = np.true_divide(image_matrix, 1000)
    image_matrix    = matts_matrix_functions.convert_matrix_to_tuple(image_matrix)
    outline         = matts_matrix_functions.get_outline_of_region(image_matrix, current_cell_outline_x_coords,current_cell_outline_y_coords)

    #Draw Position Indicator On Graph
    """
    p_trace = matts_admin_functions.create_list_of_zeros(number_of_timepoints)
    p_trace[cell_trace_timepoint] = 1000
    position_indicator_axes[current_cell_trace_plane][current_cell_trace_cell].remove()
    position_indicator_axes[current_cell_trace_plane][current_cell_trace_cell] = cell_trace_figures[current_cell_trace_plane][current_cell_trace_cell].add_subplot(111)
    position_indicator_axes[current_cell_trace_plane][current_cell_trace_cell].plot(flourescence_traces[current_cell_trace_plane][current_cell_trace_cell])
    position_indicator_axes[current_cell_trace_plane][current_cell_trace_cell].plot(np.multiply(spike_matricies[current_cell_trace_plane][current_cell_trace_cell],np.max(flourescence_traces[current_cell_trace_plane][current_cell_trace_cell])))
    position_indicator_axes[current_cell_trace_plane][current_cell_trace_cell].plot(p_trace)
    cell_trace_canvases[current_cell_trace_plane][current_cell_trace_cell].draw()
    cell_trace_canvases[current_cell_trace_plane][current_cell_trace_cell].update()
    """

    if spike_matricies[current_cell_trace_plane][current_cell_trace_cell][cell_trace_timepoint] == 1:
        colour = (0,1,0)
    else:
        colour = (1,0,0)

    image_matrix    = matts_matrix_functions.highlight_region(image_matrix , outline, colour)

    display_colour_image(main.traces_cell_figure, main.traces_cell_canvas, image_matrix)

def load_cell_trace_tiff_file(plane,file_number):
    file_name = str(folder_location + r"\\plane" + str(plane) + r"\\reg_tif\\file_chan" + str(file_number) + ".tif")
    image = io.imread(file_name)
    return image


def create_radio_function(cell):

    def radio_function():
        global current_cell_trace_file
        global current_cell_trace_cell
        global current_cell_trace_plane
        global current_cell_outline_x_coords
        global current_cell_outline_y_coords

        plane = main.traces_combox_box.currentIndex()
        current_cell_trace_plane = plane
        current_cell_trace_cell = cell
        stats_file_location = str(folder_location + str(r"\plane" + str(plane)) + "\stat.npy")
        stats = np.load(stats_file_location)

        current_cell_median = stats[cell]["med"]

        #Load Raw ROI Mask Cords
        current_cell_outline_x_coords = stats[cell]["xpix"]
        current_cell_outline_y_coords = stats[cell]["ypix"]

        #Covert These To "Window" Cords
        current_cell_outline_x_coords = current_cell_outline_x_coords - current_cell_median[1] + image_radius
        current_cell_outline_y_coords = current_cell_outline_y_coords - current_cell_median[0] + image_radius

        file_number = cell_trace_timepoint/200
        file_number = str(file_number).zfill(3)

        current_cell_trace_file = load_cell_trace_tiff_file(plane,file_number)
        update_cell_trace_tiff_index(current_cell_trace_file,neuron_coordinates[plane][cell])

        timeseries_length =  np.shape(flourescence_traces[plane][cell])
        main.traces_cell_slider.setMaximum(timeseries_length[0])



    return radio_function


def plot_flouresence_traces(plane,cell_from,cell_to):
    global cell_trace_figures
    global cell_trace_canvases

    flourescence_radiobuttons = []
    radial_functions = []

    counter = 0
    for cell in range(cell_from,cell_to):
        if neuron_decisions[plane][cell] == True:

            #Create Canvas Widget
            cell_trace_figures[plane][cell],cell_trace_canvases[plane][cell],widget = create_canvas_widget(20000,400)


            #Load and Scale Traces
            f_trace         = flourescence_traces[plane][cell]
            max_f = np.max(f_trace)
            mean_f = np.mean(f_trace)

            spike_trace     = np.multiply(spike_matricies[plane][cell], max_f)
            calcium_trace   = np.multiply(calcium_traces[plane][cell],  mean_f)
            baseline_trace  = np.multiply(baseline_traces[plane][cell], mean_f)
            fourier_trace   = np.multiply(fourier_traces[plane][cell], mean_f)

            #Plot Traces On Canvas
            plot_trace_graph([spike_trace,f_trace,baseline_trace,calcium_trace,fourier_trace],cell_trace_figures[plane][cell],cell_trace_canvases[plane][cell],"Neuron: " + str(cell),plane,cell)

            #Create Radio Buttons
            flourescence_radiobuttons.append(QRadioButton(str(cell)))
            radio_function = create_radio_function(cell)
            radial_functions.append(radio_function)
            flourescence_radiobuttons[counter].toggled.connect(radio_function)

            main.traces_scroll_area_contents_layout.addWidget(widget,                                   counter, 1)
            main.traces_scroll_area_contents_layout.addWidget(flourescence_radiobuttons[counter],       counter, 0)

            counter += 1
            app.processEvents()



def plot_line_graph(y_cords_list,figure,canvas,title):
    figure.clear()
    axis = figure.add_subplot(111)
    axis.set_title(str(title))
    #axis.set_ylim(0,1.2)
    #axis.set_xlim(-1.5, 1.5

    axis.plot(y_cords_list)

    canvas.draw()
    canvas.update()


def plot_multiple_lines(traces,figure,canvas,title):
    figure.clear()
    axis = figure.add_subplot(111)
    axis.set_title(str(title))

    for trace in traces:
        axis.plot(trace)

    canvas.draw()
    canvas.update()


def plot_trace_graph(traces,figure,canvas,title,plane,cell):
    global position_indicator_axes

    #xticks(np.arange(0, len(traces[0]), step=100))

    figure.clear()
    figure.frameon = False
    axis = figure.add_subplot(111)
    axis.set_title(str(title))
    axis.set_xticks(np.arange(0, len(traces[0]), step=250))
    axis.tick_params(labelsize=6)


    colours = ['ro',"c","k","m","g"]

    for trace_index in range(len(traces)):
        axis.plot(traces[trace_index],colours[trace_index])

    #axis.plot(traces[0][0],traces[0][1],colours[0])

    #position_indicator_axes[plane][cell] = figure.add_subplot(111)

    canvas.draw()
    canvas.update()


def plot_histogram(figure,canvas,list,title):

    figure.clear()
    axis = figure.add_subplot(111)
    histogram, bin_edges = np.histogram(list)
    n, bins, patches = axis.hist(x=list, bins='auto', color='#0504aa', alpha=1, rwidth=0.85)
    axis.grid(axis='y', alpha=0.75)
    axis.set_ylabel('Frequency')
    axis.set_title(title)

    axis.set_xticks(np.arange(0, np.max(list), step=200))
    canvas.draw()
    canvas.update()



class Window(QDialog):

    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle("Calcium Imaging Pipeline")
        #self.setGeometry(100,100,1500,700)
        self.showMaximized()

        self.create_main_layout()
        self.create_region_selection_tab()
        self.create_traces_tab()
        self.create_minimum_value_histogram_tab()

    def create_main_layout(self):
        self.main_layout = QGridLayout()
        self.setLayout(self.main_layout)

        self.tabs = QTabWidget()
        self.main_layout.addWidget(self.tabs)

    def threshold_slider_value_changed(self):
        global cell_probability_threshold
        cell_probability_threshold = float(main.probabiltiy_slider.value())/100
        main.probabiltiy_slider_value_label.setText(str(cell_probability_threshold))

    def show_traces_button_clicked(self):

        for i in reversed(range(main.traces_scroll_area_contents_layout.count())):
            main.traces_scroll_area_contents_layout.itemAt(i).widget().setParent(None)

        selected_plane = main.traces_combox_box.currentIndex()
        cell_from = main.neuron_start_spinner.value()
        cell_to = main.neuron_end_spinner.value()
        plot_flouresence_traces(selected_plane,cell_from,cell_to)

    def show_all_traces_button_clicked(self):

        for i in reversed(range(main.traces_scroll_area_contents_layout.count())):
            main.traces_scroll_area_contents_layout.itemAt(i).widget().setParent(None)

        selected_plane = main.traces_combox_box.currentIndex()
        cell_from = 0
        cell_to = len(neuron_decisions[selected_plane])

        plot_flouresence_traces(selected_plane, cell_from, cell_to)


    def trace_combo_box_changed(self):
        selected_plane = main.traces_combox_box.currentIndex()
        max = len(flourescence_traces[selected_plane])

        main.neuron_end_spinner.setMaximum(max)
        main.neuron_start_spinner.setMaximum(max)

    def threshold_minimum_value_clicked(self):
        global cell_noise_threshold
        cell_noise_threshold = main.minimum_value_input.value()
        threshold_cell_noises()

    def set_number_of_planes(self):
        global number_of_planes
        number_of_planes = self.number_of_planes_spinner.value()

    def cell_trace_slider_changed(self):
        global cell_trace_timepoint
        global current_cell_trace_file
        global current_cell_trace_cell

        plane = main.traces_combox_box.currentIndex()

        current_file_number = cell_trace_timepoint / 200

        cell_trace_timepoint = self.traces_cell_slider.value()
        self.traces_cell_slider_label.setText(str(cell_trace_timepoint))

        new_file_number = cell_trace_timepoint / 200

        coordinates = neuron_coordinates[plane][current_cell_trace_cell]

        #If we're still in the same file, just update the index
        if new_file_number == current_file_number:
            update_cell_trace_tiff_index(current_cell_trace_file,coordinates)

        else:
            new_file_number = str(new_file_number).zfill(3)
            current_cell_trace_file = load_cell_trace_tiff_file(plane, new_file_number)
            update_cell_trace_tiff_index(current_cell_trace_file,coordinates)

    def create_region_combo_box_function(self,plane):

        def combo_box_function():
            global currently_selected_region
            region = str(region_combo_boxes[plane].currentText())
            currently_selected_region[plane] = region

        return combo_box_function

    def parameters_changed(self):
        global wdt
        global lamb
        global varB
        global varC
        global initial_baseline
        global Cmean
        global fourier_cutoff

        wdt                 = float(self.wdt_slider.value())                 / 100
        lamb                = float(self.lamb_slider.value())                / 100
        varB                = float(self.varB_slider.value())                / 10000
        varC                = float(self.varC_slider.value())                / 100
        initial_baseline    = float(self.initial_baseline_slider.value())    / 100
        Cmean               = float(self.Cmean_slider.value())               / 100
        fourier_cutoff      = self.fourier_cutoff_slider.value()

        self.wdt_slider_value_label                 .setText(str(wdt))
        self.lamb_slider_value_label                .setText(str(lamb))
        self.varB_slider_value_label                .setText(str(varB))
        self.varC_slider_value_label                .setText(str(varC))
        self.initial_baseline_slider_value_label    .setText(str(initial_baseline))
        self.Cmean_slider_value_label               .setText(str(Cmean))
        self.fourier_cutoff_slider_value_label      .setText(str(fourier_cutoff))


    def populate_region_selection_tab(self):
        self.region_figures = matts_admin_functions.create_empty_list(number_of_planes)
        self.region_canvases = matts_admin_functions.create_empty_list(number_of_planes)
        self.region_widgets = matts_admin_functions.create_empty_list(number_of_planes)
        self.roi_displays = matts_admin_functions.create_empty_list(number_of_planes)
        global region_combo_boxes

        for plane in range(number_of_planes):
            figure, canvas, widget = create_canvas_widget(400, 400)
            self.region_figures[plane] = figure
            self.region_canvases[plane] = canvas
            self.region_widgets[plane] = widget

            self.roi_displays[plane] = QGraphicsView()
            self.roi_displays[plane].setMinimumWidth(500)
            region_combo_boxes[plane] = QComboBox()

            region_combo_boxes[plane].addItem("Forebrain")
            region_combo_boxes[plane].addItem("Tectum")
            region_combo_boxes[plane].addItem("Midbrain")
            region_combo_boxes[plane].addItem("Hindbrain")

            combo_box_function = self.create_region_combo_box_function(plane)
            region_combo_boxes[plane].currentIndexChanged.connect(combo_box_function)

            self.region_selection_scroll_area_contents_layout.addWidget(region_combo_boxes[plane], plane, 0)
            self.region_selection_scroll_area_contents_layout.addWidget(self.region_widgets[plane], plane, 1)
            self.region_selection_scroll_area_contents_layout.addWidget(self.roi_displays[plane], plane, 2)


    def create_region_selection_tab(self):

        #Create Tab and Add It To Tab Widet
        self.region_selection_tab = QWidget()
        self.region_selection_tab_layout = QGridLayout()
        self.region_selection_tab.setLayout(self.region_selection_tab_layout)

        #Input Directory Button
        self.input_directory_button = QPushButton("Load Suite2p Data")
        self.input_directory_button.clicked.connect(load_data)

        #Number of Planes Input
        self.number_of_planes_spinner = QSpinBox()
        self.number_of_planes_spinner.setValue(number_of_planes)
        self.number_of_planes_label = QLabel("Number Of Planes: ")
        self.number_of_planes_spinner.valueChanged.connect(self.set_number_of_planes)

        #Cell Probability Q-Slider
        self.probabiltiy_slider = QSlider(Qt.Horizontal)
        self.probabiltiy_slider.valueChanged.connect(self.threshold_slider_value_changed)
        self.probabiltiy_slider.setValue(cell_probability_threshold*100)
        self.probability_threshold_button = QPushButton("Threshold Cell Probability")
        self.probability_threshold_button.clicked.connect(threshold_cell_probability)
        self.probabiltiy_slider_label = QLabel("Cell Probability Threshold")
        self.probabiltiy_slider_value_label = QLabel(str(cell_probability_threshold))

        #Save Traces Button
        self.save_traces_button = QPushButton("Save Traces")
        self.save_traces_button.clicked.connect(save_traces)

        #Extract Spikes Button
        self.extract_spikes_button = QPushButton("Extract Spikes")
        self.extract_spikes_button.clicked.connect(extract_spikes)

        #Number of Neurons Button and Label
        self.number_of_neurons_button = QPushButton("Calculate Number of Remaining Neurons")
        self.number_of_neurons_button.clicked.connect(calculate_number_of_cells)
        self.number_of_neurons_label = QLabel()


        #BCL Parameter Tuning
        self.bcl_parameter_label = QLabel("BCL Parameters")

        self.wdt_slider                 = QSlider(Qt.Horizontal)
        self.lamb_slider                = QSlider(Qt.Horizontal)
        self.varB_slider                = QSlider(Qt.Horizontal)
        self.varC_slider                = QSlider(Qt.Horizontal)
        self.initial_baseline_slider    = QSlider(Qt.Horizontal)
        self.Cmean_slider               = QSlider(Qt.Horizontal)
        self.fourier_cutoff_slider      = QSlider(Qt.Horizontal)

        self.wdt_slider_label               = QLabel("Spiking Probability")
        self.lamb_slider_label              = QLabel("Calcium Decay")
        self.varB_slider_label              = QLabel("Baseline Variance")
        self.varC_slider_label              = QLabel("Calcium Variance")
        self.initial_baseline_slider_label  = QLabel("Initial Baseline")
        self.Cmean_slider_label             = QLabel("Mean Calcium")
        self.fourier_cutoff_slider_label    = QLabel("Fourier Cutoff")

        self.wdt_slider                 .setTickInterval(0.1    * 100)
        self.lamb_slider                .setTickInterval(0.1    * 100)
        self.varB_slider                .setTickInterval(0.0001 * 10000)
        self.varC_slider                .setTickInterval(0.1    * 100)
        self.initial_baseline_slider    .setTickInterval(0.1    * 100)
        self.Cmean_slider               .setTickInterval(0.1    * 100)
        self.fourier_cutoff_slider      .setTickInterval(100)

        self.wdt_slider                 .setMaximum(1   * 100)
        self.lamb_slider                .setMaximum(10  * 100)
        self.varB_slider                .setMaximum(1   * 10000)
        self.varC_slider                .setMaximum(5   * 100)
        self.initial_baseline_slider    .setMaximum(5   * 100)
        self.Cmean_slider               .setMaximum(5   * 100)
        self.fourier_cutoff_slider      .setMaximum(30000)

        self.wdt_slider                 .setValue(wdt   * 100)
        self.lamb_slider                .setValue(lamb  * 100)
        self.varB_slider                .setValue(varB  * 10000)
        self.varC_slider                .setValue(varC  * 100)
        self.initial_baseline_slider    .setValue(initial_baseline * 100)
        self.Cmean_slider               .setValue(Cmean * 100)
        self.fourier_cutoff_slider      .setValue(fourier_cutoff)

        self.wdt_slider_value_label               = QLabel(str(wdt))
        self.lamb_slider_value_label              = QLabel(str(lamb))
        self.varB_slider_value_label              = QLabel(str(varB))
        self.varC_slider_value_label              = QLabel(str(varC))
        self.initial_baseline_slider_value_label  = QLabel(str(initial_baseline))
        self.Cmean_slider_value_label             = QLabel(str(Cmean))
        self.fourier_cutoff_slider_value_label    = QLabel(str(fourier_cutoff))

        self.wdt_slider                 .valueChanged.connect(self.parameters_changed)
        self.lamb_slider                .valueChanged.connect(self.parameters_changed)
        self.varB_slider                .valueChanged.connect(self.parameters_changed)
        self.varC_slider                .valueChanged.connect(self.parameters_changed)
        self.initial_baseline_slider    .valueChanged.connect(self.parameters_changed)
        self.Cmean_slider               .valueChanged.connect(self.parameters_changed)
        self.fourier_cutoff_slider      .valueChanged.connect(self.parameters_changed)




        #Create Scroll Area
        self.region_selection_scroll_area_contents = QWidget()
        self.region_selection_scroll_area_contents_layout = QGridLayout()
        self.region_selection_scroll_area_contents.setLayout(self.region_selection_scroll_area_contents_layout)
        self.region_selection_scroll_area = QScrollArea()
        self.region_selection_scroll_area.setWidget(self.region_selection_scroll_area_contents)
        self.region_selection_scroll_area.setWidgetResizable(True)

        #Add Widgets To Layout
        self.region_selection_tab_layout.setColumnMinimumWidth(0, 1500)
        self.region_selection_tab_layout.addWidget(self.input_directory_button,         0, 1, 01, 3)
        self.region_selection_tab_layout.addWidget(self.number_of_planes_spinner,       1, 3, 01, 1)
        self.region_selection_tab_layout.addWidget(self.number_of_planes_label,         1, 1, 01, 2)
        self.region_selection_tab_layout.addWidget(self.region_selection_scroll_area,   0, 0, 20, 1)
        self.region_selection_tab_layout.addWidget(self.probabiltiy_slider_label,       2, 1, 01, 1)
        self.region_selection_tab_layout.addWidget(self.probabiltiy_slider,             2, 2, 01, 1)
        self.region_selection_tab_layout.addWidget(self.probabiltiy_slider_value_label, 2, 3, 01, 1)
        self.region_selection_tab_layout.addWidget(self.probability_threshold_button,   3, 1, 01, 3)
        self.region_selection_tab_layout.addWidget(self.save_traces_button,             4, 1, 01, 3)
        self.region_selection_tab_layout.addWidget(self.extract_spikes_button,          5, 1, 01, 3)
        self.region_selection_tab_layout.addWidget(self.number_of_neurons_button,       6, 1, 01, 1)
        self.region_selection_tab_layout.addWidget(self.number_of_neurons_label,        6, 2, 01, 2)

        self.region_selection_tab_layout.addWidget(self.bcl_parameter_label,                     7, 01, 01, 03)

        self.region_selection_tab_layout.addWidget(self.wdt_slider_label,                        8, 01, 01, 01)
        self.region_selection_tab_layout.addWidget(self.lamb_slider_label,                       9, 01, 01, 01)
        self.region_selection_tab_layout.addWidget(self.varB_slider_label,                      10, 01, 01, 01)
        self.region_selection_tab_layout.addWidget(self.varC_slider_label,                      11, 01, 01, 01)
        self.region_selection_tab_layout.addWidget(self.initial_baseline_slider_label,          12, 01, 01, 01)
        self.region_selection_tab_layout.addWidget(self.Cmean_slider_label,                     13, 01, 01, 01)
        self.region_selection_tab_layout.addWidget(self.fourier_cutoff_slider_label,            14, 01, 01, 01)

        self.region_selection_tab_layout.addWidget(self.wdt_slider,                              8, 02, 01, 01)
        self.region_selection_tab_layout.addWidget(self.lamb_slider,                             9, 02, 01, 01)
        self.region_selection_tab_layout.addWidget(self.varB_slider,                            10, 02, 01, 01)
        self.region_selection_tab_layout.addWidget(self.varC_slider ,                           11, 02, 01, 01)
        self.region_selection_tab_layout.addWidget(self.initial_baseline_slider,                12, 02, 01, 01)
        self.region_selection_tab_layout.addWidget(self.Cmean_slider,                           13, 02, 01, 01)
        self.region_selection_tab_layout.addWidget(self.fourier_cutoff_slider,                  14, 02, 01, 01)

        self.region_selection_tab_layout.addWidget(self.wdt_slider_value_label,                  8, 03, 01, 01)
        self.region_selection_tab_layout.addWidget(self.lamb_slider_value_label,                 9, 03, 01, 01)
        self.region_selection_tab_layout.addWidget(self.varB_slider_value_label,                10, 03, 01, 01)
        self.region_selection_tab_layout.addWidget(self.varC_slider_value_label,                11, 03, 01, 01)
        self.region_selection_tab_layout.addWidget(self.initial_baseline_slider_value_label,    12, 03, 01, 01)
        self.region_selection_tab_layout.addWidget(self.Cmean_slider_value_label,               13, 03, 01, 01)
        self.region_selection_tab_layout.addWidget(self.fourier_cutoff_slider_value_label,      14, 03, 01, 01)



        self.tabs.addTab(self.region_selection_tab, "Region Selection")


    def create_traces_tab(self):
        self.traces_tab = QWidget()
        self.traces_tab_layout = QGridLayout()
        self.traces_tab.setLayout(self.traces_tab_layout)

        self.traces_cell_figure, self.traces_cell_canvas, self.traces_cell_widget = create_canvas_widget(400,400)
        self.traces_cell_canvas_label = QLabel("Registered Tiff")

        self.traces_cell_slider = QSlider(Qt.Horizontal)
        self.traces_cell_slider.valueChanged.connect(self.cell_trace_slider_changed)
        self.traces_cell_slider_label = QLabel("0")
        self.traces_cell_slider.setFixedWidth(400)

        self.traces_scroll_area_contents = QWidget()
        self.traces_scroll_area_contents_layout = QGridLayout()
        self.traces_scroll_area_contents.setLayout(self.traces_scroll_area_contents_layout)

        self.traces_scroll_area = QScrollArea()
        self.traces_scroll_area.setWidget(self.traces_scroll_area_contents)
        self.traces_scroll_area.setWidgetResizable(True)

        self.traces_combox_box = QComboBox()
        self.traces_combox_box.currentIndexChanged.connect(self.trace_combo_box_changed)

        self.neuron_start_spinner = QSpinBox()
        self.neuron_end_spinner = QSpinBox()
        self.neuron_end_spinner.setValue(10)
        self.neuron_start_spinner_label = QLabel("From Neuron:")
        self.neuron_end_spinner_label = QLabel("To Neuron:")

        self.show_traces_button = QPushButton("Display Traces")
        self.show_traces_button.clicked.connect(self.show_traces_button_clicked)

        self.show_all_traces_button = QPushButton("Display All Traces ")
        self.show_all_traces_button.clicked.connect(self.show_all_traces_button_clicked)

        self.traces_tab_layout.addWidget(self.traces_cell_canvas_label,     0, 0, 1, 1)
        self.traces_tab_layout.addWidget(self.traces_cell_widget,           1, 0, 1, 1)
        self.traces_tab_layout.addWidget(self.traces_cell_slider,        2, 0, 1, 1)
        self.traces_tab_layout.addWidget(self.traces_cell_slider_label,  3, 0, 1, 1)
        self.traces_tab_layout.addWidget(self.traces_scroll_area,           0, 1, 20, 1)
        self.traces_tab_layout.addWidget(self.traces_combox_box,            0, 2, 1, 2)

        self.traces_tab_layout.addWidget(self.neuron_start_spinner,         1, 3, 1, 1)
        self.traces_tab_layout.addWidget(self.neuron_start_spinner_label,   1, 2, 1, 1)

        self.traces_tab_layout.addWidget(self.neuron_end_spinner,           2, 3, 1, 1)
        self.traces_tab_layout.addWidget(self.neuron_end_spinner_label,     2, 2, 1, 1)

        self.traces_tab_layout.addWidget(self.show_traces_button,           3, 2, 1, 2)
        self.traces_tab_layout.addWidget(self.show_all_traces_button,       4, 2, 1, 2)

        self.tabs.addTab(self.traces_tab, "Traces")


    def create_minimum_value_histogram_tab(self):
        self.minimum_values_tab = QWidget()
        self.minimum_values_tab_layout = QGridLayout()
        self.minimum_values_tab.setLayout(self.minimum_values_tab_layout)

        self.minimum_values_figure,self.minimum_values_canvas,self.minimum_values_widget = create_canvas_widget(800,800)

        self.calculate_minimum_values_button = QPushButton("Calculate Highest Minimums")
        self.calculate_minimum_values_button.clicked.connect(calculate_highest_minimums)

        self.minimum_values_progress_bar = QProgressBar()

        self.minimum_value_input = QSpinBox()
        self.minimum_value_input.setMaximum(10000)
        self.minimum_value_input.setValue(50)
        self.minimum_value_input_label = QLabel("Minimum Value Threshold")


        self.minimum_value_threshold_button = QPushButton("Threshold Minimum Value")
        self.minimum_value_threshold_button.clicked.connect(self.threshold_minimum_value_clicked)

        self.minimum_values_tab_layout.addWidget(self.minimum_values_canvas,            0, 0, 1, 3)
        self.minimum_values_tab_layout.addWidget(self.calculate_minimum_values_button,  1, 0, 1, 3)
        self.minimum_values_tab_layout.addWidget(self.minimum_values_progress_bar,      2, 0, 1, 3)
        self.minimum_values_tab_layout.addWidget(self.minimum_value_threshold_button,   3, 0, 1, 1)
        self.minimum_values_tab_layout.addWidget(self.minimum_value_input_label,        3, 1, 1, 1)
        self.minimum_values_tab_layout.addWidget(self.minimum_value_input,              3, 2, 1, 1)
        self.tabs.addTab(self.minimum_values_tab,"Highest Minimum Value Histogram")


def save_traces():

    output_file_folder = str(QFileDialog.getExistingDirectory(main, "Open a folder"))

    #Split Traces and ROIs by Region
    forebrain_traces        = []
    forebrain_coordinates   = []
    tectum_traces           = []
    tectum_coordinates      = []
    midbrain_traces         = []
    midbrain_coordinates    = []
    hindbrain_traces        = []
    hindbrain_coordinates   = []

    counter = 0
    for plane in range(number_of_planes):
        for cell in range(number_of_neurons[plane]):
            if neuron_decisions[plane][cell] == True:

                location = neuron_locations[plane][cell]

                if location == "Forebrain":
                    forebrain_traces.append(spike_matricies[plane][cell])
                    forebrain_coordinates.append(neuron_coordinates[plane][cell])

                elif location == "Tectum":
                    tectum_traces.append(spike_matricies[plane][cell])
                    tectum_coordinates.append(neuron_coordinates[plane][cell])

                elif location == "Midbrain":
                    midbrain_traces.append(spike_matricies[plane][cell])
                    midbrain_coordinates.append(neuron_coordinates[plane][cell])

                elif location == "Hindbrain":
                    hindbrain_traces.append(spike_matricies[plane][cell])
                    hindbrain_coordinates.append(neuron_coordinates[plane][cell])

                counter += 1

    #Turn Into Numpy Arrays
    forebrain_traces        = np.array(forebrain_traces)
    forebrain_coordinates   = np.array(forebrain_coordinates)
    tectum_traces           = np.array(tectum_traces)
    tectum_coordinates      = np.array(tectum_coordinates)
    midbrain_traces         = np.array(midbrain_traces)
    midbrain_coordinates    = np.array(midbrain_coordinates)
    hindbrain_traces        = np.array(hindbrain_traces)
    hindbrain_coordinates   = np.array(hindbrain_coordinates)

    #Save Binary Matricies
    np.savetxt(output_file_folder + "\\Forebrain Binary Matrix.csv",  forebrain_traces,       delimiter=",")
    np.savetxt(output_file_folder + "\\Tectum Binary Matrix.csv",     tectum_traces,          delimiter=",")
    np.savetxt(output_file_folder + "\\Midbrain Binary Matrix.csv",   midbrain_traces,        delimiter=",")
    np.savetxt(output_file_folder + "\\Hindbrain Binary Matrix.csv",  hindbrain_traces,       delimiter=",")

    #Save ROI Coordinates
    np.savetxt(output_file_folder + "\\Forebrain ROIs.csv",           forebrain_coordinates,  delimiter=",")
    np.savetxt(output_file_folder + "\\Tectum ROIs.csv",              tectum_coordinates,     delimiter=",")
    np.savetxt(output_file_folder + "\\Midbrain ROIs.csv",            midbrain_coordinates,   delimiter=",")
    np.savetxt(output_file_folder + "\\Hindbrain ROIs.csv",           hindbrain_coordinates,  delimiter=",")

    print "saved!"



def threshold_cell_probability():
    for plane in range(number_of_planes):
        for cell in range(len(cell_probabilities[plane])):
            if cell_probabilities[plane][cell] > cell_probability_threshold:
                cell_probability_decisions[plane][cell] = True
            else:
                cell_probability_decisions[plane][cell] = False

    threshold_cells()
    draw_all_roi_scenes()

def threshold_cell_noises():
    for plane in range(number_of_planes):
        for cell in range(len(cell_noises[plane])):

            if cell_noises[plane][cell] > cell_noise_threshold:
                cell_noise_decisions[plane][cell] = True
            else:
                cell_noise_decisions[plane][cell] = False


    threshold_cells()
    draw_all_roi_scenes()

def create_canvas_widget(width,height):
    figure = Figure()
    canvas = FigureCanvas(figure)
    widget = QWidget()
    layout = QGridLayout()
    widget.setLayout(layout)
    layout.addWidget(canvas)
    widget.setFixedSize(width, height)
    return figure, canvas, widget


def load_mean_images():
    global region_images
    global tiff_dimensions

    region_images = []
    for plane in range(number_of_planes):
        ops_file_location = str(folder_location + str(r"\plane" + str(plane)) + "\ops.npy")
        ops = np.load(ops_file_location)
        ops = ops[()]
        region_images.append(ops["meanImg"])
        display_image(main.region_figures[plane],main.region_canvases[plane],region_images[plane])

    tiff_dimensions = [region_images[0].shape[1], region_images[0].shape[0]]


def display_image(figure,canvas,image):
    figure.clear()
    axis = figure.add_subplot(111)
    axis.imshow(image, cmap='gray')
    axis.axis("off")
    plt.show()
    canvas.draw()
    canvas.update()

def display_colour_image(figure,canvas,image):
    figure.clear()
    axis = figure.add_subplot(111)
    axis.imshow(image, vmin=0,vmax=1, cmap="gray")
    axis.axis("off")
    plt.show()
    canvas.draw()
    canvas.update()


def create_pseduo_rois(number):
    global neuron_decisions

    coordinates = []

    for roi in range(number):
        x = np.random.uniform(0,tiff_dimensions[0])
        y = np.random.uniform(0,tiff_dimensions[1])
        coordinates.append((x,y))

    neuron_decisions = matts_admin_functions.create_list_of_item(number,True)

    return coordinates


def draw_all_roi_scenes():
    for plane in range(number_of_planes):
        draw_roi_scene(plane)


def draw_roi_scene(plane):
    scene = QGraphicsScene()
    scale_factor = 3
    #scene.setSceneRect(0,0, tiff_dimensions[0], tiff_dimensions[1])

    roi_size = 10
    roi_elipses = matts_admin_functions.create_empty_list(len(neuron_coordinates[plane]))

    y_offset = -50
    x_offset = -50

    for roi_index in range(len(neuron_coordinates[plane])):

        x_cord = neuron_coordinates[plane][roi_index][1]
        y_cord = neuron_coordinates[plane][roi_index][0]


        probability = cell_probabilities[plane][roi_index]

        #Determine Colour
        if neuron_locations[plane][roi_index] == 0:
            colour = get_colour(probability,"hot",0.8)

        elif neuron_locations[plane][roi_index] == "Forebrain":
            colour = get_colour(probability, "Blues", 0.8)

        elif neuron_locations[plane][roi_index] == "Tectum":
            colour = get_colour(probability, "Purples", 0.8)

        elif neuron_locations[plane][roi_index] == "Midbrain":
            colour = get_colour(probability, "Greens", 0.8)

        elif neuron_locations[plane][roi_index] == "Hindbrain":
            colour = get_colour(probability, "Reds", 0.8)


        roi_elipses[roi_index] = QGraphicsEllipseItem((x_cord + x_offset), (y_cord + y_offset), roi_size, roi_size)
        roi_elipses[roi_index].setBrush(QBrush(colour))
        roi_elipses[roi_index].setPos(x_cord / scale_factor, y_cord / scale_factor)

        if neuron_decisions[plane][roi_index] == False:
            roi_elipses[roi_index].setOpacity(0.2)

        scene.addItem(roi_elipses[roi_index])


    main.roi_displays[plane].setScene(scene)



def get_colour(input_value,colour_map,scale_factor):
    input_value = input_value * scale_factor
    cmap = cm.get_cmap(colour_map)
    float_tuple = cmap(input_value)
    matplot_to_q_colour_conversion_factor = 255
    colour = QColor(float_tuple[0]*matplot_to_q_colour_conversion_factor,
                          float_tuple[1]*matplot_to_q_colour_conversion_factor,
                          float_tuple[2]*matplot_to_q_colour_conversion_factor)

    return colour



def create_lasso_functions():

    plane_functions = matts_admin_functions.create_empty_list(number_of_planes)

    for plane in range(number_of_planes):
        plane_functions[plane] = create_plane_function(plane)

    return plane_functions



def create_plane_function(plane):

    def plane_selected(verts):
        global currently_selected_region

        x_tiff, y_tiff = np.meshgrid(np.arange(region_images[plane].shape[1]), np.arange(region_images[plane].shape[0]))
        pix = np.vstack((x_tiff.flatten(), y_tiff.flatten())).T
        p = Path(verts)
        ind = p.contains_points(pix, radius=1)
        selected = np.zeros(tiff_dimensions)
        selected.flat[ind] = 1

        for neuron in range(len(neuron_decisions[plane])):
            y_cord = int(neuron_coordinates[plane][neuron][0])
            x_cord = int(neuron_coordinates[plane][neuron][1])

            if selected[y_cord][x_cord] == 1:
                within_region_decisions[plane][neuron] = True
                neuron_locations[plane][neuron] = currently_selected_region[plane]
            #else:
                #within_region_decisions[plane][neuron] = False



        threshold_cells()
        draw_roi_scene(plane)
        app.processEvents()

    return plane_selected


def create_lassos():
    lassos = matts_admin_functions.create_empty_list(number_of_planes)
    lasso_functions = create_lasso_functions()

    for plane in range(number_of_planes):
        axis = main.region_figures[plane].get_axes()
        lassos[plane] = LassoSelector(axis[0], lasso_functions[plane])

    return lassos

def extract_spikes():
    global spike_matricies
    global calcium_traces
    global baseline_traces

    counter = 0
    for plane in range(number_of_planes):
        for neuron in range(len(neuron_decisions[plane])):
            if neuron_decisions[plane][neuron] == True:
                print counter
                average_f_trace = average_trace(flourescence_traces[plane][neuron],4)
                calcium_signal, baseline_signal, spikes, fourier = PCL_Parameters.bcl_function_parameters(average_f_trace,wdt,lamb,varB,varC,initial_baseline,Cmean,frequency,fourier_cutoff)
                spike_matricies [plane][neuron] = spikes
                calcium_traces  [plane][neuron] = calcium_signal
                baseline_traces [plane][neuron] = baseline_signal
                fourier_traces  [plane][neuron] = fourier
                counter += 1

if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    main.show()
    sys.exit(app.exec_())
