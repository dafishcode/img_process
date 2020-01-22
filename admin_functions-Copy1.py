def name(path, experiment, f): #add 0 to number if above 10
    os.chdir(path + 'Project/' + experiment)
    if f < 10:
        num = '0' + str(f)
    else:
        num = str(f)
    return(num)