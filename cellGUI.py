import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from matplotlib.figure import Figure
from PIL import Image, ImageSequence
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas



class Window(QWidget):

    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle("To cell or not to cell")
        self.setGeometry(100, 100, 1500, 700)

        self.create_main_layout()

    def create_main_layout(self):
        self.layout = QGridLayout()
        self.setLayout(self.layout)
        self.cellcombobox = QComboBox()
        self.slider = QSlider(Qt.Horizontal)
        self.button = QPushButton('pushme')
        self.layout.addWidget(self.slider, 1, 1)
        self.layout.addWidget(self.button, 5, 3)
        self.layout.addWidget(self.cellcombobox, 5, 4)
        self.image = Figure()
        self.imagecanvas = FigureCanvas(self.image)
        self.layout.addWidget(self.imagecanvas)



    #make function - add cells
    for cell in range leng etccc
        self.cellcombobox.addItem(string - cell (ie number of 0 to len noise cll))

    function - when selected on combo box load tiff

        self.cellcombobox.currentindexchanged.connect(function name)

    make function called when combobox is changed -



    global variable - currently select3d cel\

    what cell we want
    self.combobox.currentindex (noise cell £number[2] ie choose plane)

    load tiffs - as reg tiffs
    io.imread - etc


    update figure canvas

    def display_image(figure, canvas, image):
        figure.clear()
        axis = figure.add_subplot(111)
        axis.imshow(image, cmap='gray')
        axis.axis("off")
        plt.show()
        canvas.draw()
        canvas.update()
        
    def display image - from formatter


    find exactly where cell is in image -

    whenever slider value changes, function to change image []value


    global variable - which cell is currently collected (exists in between functions)


    button - current index from current box

    add to kepttr


    # tempalte to make window


app = QApplication(sys.argv)

thing = Window()
thing.show()

sys.exit(app.exec_())
# tiffdisplayfigure = Figure()
# tiffdisplaycanvas = FigureCanvas(tiffdisplayfigure)


# make function

# tiffdisplayfigure.clear()
# axis = figure.add_subplot(111)
# axis.imshow(practice[0], cmap = 'grey')
# plt.show()
# tiffdisplaycanvas.draw()
# tiffdisplaycanvas.update()