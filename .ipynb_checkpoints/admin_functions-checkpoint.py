from matplotlib import cm
import random
import matplotlib.pyplot as plt


def split_list_by_spaces(list):
    number_list = []
    list_string = str(list)
    split_list_string = list_string.split(" ")

    for item in split_list_string:
        item = item.replace("[", "")
        item = item.replace("]", "")

        if item == "":
            pass

        else:
            number = float(item)
            number_list.append(number)

    return number_list

def plot_graph(values):
    plt.plot(values)
    plt.show()



def get_boolean_with_probability(probability):
    number_1 = random.uniform(0,1)

    if probability == 0:
        return False

    if probability == 1:
        return True

    if probability > number_1:
        return True
    else:
        return False



def create_empty_array(size):
    new_array = create_empty_list(size)

    for row in range(size):
        columns = create_empty_list(size)
        new_array[row] = columns

    return new_array

def create_empty_list(size):
    new_list = []
    for number in range(size):
        new_list.append([])

    return new_list

def create_list_of_zeros(size):
    new_list = []
    for number in range(size):
        new_list.append(0)

    return new_list


def create_list_of_item(size,item):
    new_list = []
    for number in range(size):
        new_list.append(item)

    return new_list

def assign_colour(item,number_of_items,colour_map):
    colour_map = cm.get_cmap(colour_map)
    value = float(item)/number_of_items
    float_tuple = colour_map(value)
    return float_tuple



