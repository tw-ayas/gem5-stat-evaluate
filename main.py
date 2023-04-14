import math
import sys
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as tck
from benchmarks import list_benchmarks, size, depth, excel_labels

complete_dict = {}
requestor = {
    '5': 'cpu0.inst',
    '6': 'cpu0.data',
    '7': 'cpu0.dcache.prefetcher',
    '11': 'cpu1.inst',
    '12': 'cpu1.data',
    '13': 'cpu1.dcache.prefetcher',
    'total': 'total',
}


def stat_name_to_list(names, count, index):
    listD = {}
    if len(names) == index:
        return count
    indexF = index + 1
    listD[names[index]] = stat_name_to_list(names, count, indexF)
    return listD


def dict_merge(dct, merge_dct):
    for k, v in merge_dct.items():
        if k in dct and isinstance(dct[k], dict) and isinstance(merge_dct[k], dict):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]


def is_in(word, list):
    for x in list:
        if word.lower() in x.lower():
            return True
    return False


def compare_count_name(word):
    return word.lower().replace("countmin", "").replace("_", "")


def find_similar_count_min(key, data, data_orig, name_list):
    if isinstance(data, str):
        # print(name_list)
        # print(f"{key} => {data}")
        return
    for k, v in data.items():
        name_list_orig = name_list[:]
        if isinstance(data[k], dict) and "::" not in k:
            name_list.append(k)
            find_similar_count_min(k, v, data_orig, name_list)
            name_list = name_list_orig[:]
        else:

            if "countmin" in k.lower() or is_in("countmin", name_list):
                if "::" in k:
                    try:
                        key_pre_pos = k.split("::")
                        if len(key_pre_pos) > 1 and key_pre_pos[1].isnumeric():
                            k_v = requestor[key_pre_pos[1]].split(".")
                            k = key_pre_pos[0] + "::" + k_v[0]
                            if len(k_v) > 1:
                                v = {k_v[1]: v}
                                if len(k_v) > 2:
                                    v = {k_v[1]: {k_v[2]: v[k_v[1]]}}

                    except ValueError:
                        k = k
                    except KeyError:
                        k = k
                name = ".".join(name_list) + "." + k
                to_compare = data_orig.copy()
                for n in name_list[1:]:
                    if "countMin" in n:
                        n = n.replace("countMin", "").replace("_", "")
                    to_compare = to_compare[n]
                for key, value in to_compare.items():
                    name_compare = compare_count_name(k)
                    if name_compare == key.replace("_", "").lower():
                        if isinstance(v, dict):
                            for key_parent, value_parent in v.items():
                                name_count = name + "." + key_parent
                                value_parsed = value_parent
                                if isinstance(value_parsed, dict):
                                    name_count = name_count + "." + list(value_parent.items())[0][0]
                                    value_parsed = list(value_parent.items())[0][1]
                                name_count = name_count.replace('cpu0', 'cpu')
                                if name_count not in complete_dict:
                                    complete_dict[name_count] = []
                                    complete_dict[name_count].append(name_count)
                                complete_dict[name_count].append(value_parsed)
                                # print(name_count + " => " + value_parsed)
                                for key_child, value_child in value.items():
                                    compare_name = compare_count_name(k + "." + key_parent)
                                    # print("                 compare_name: " + compare_name)
                                    # print("                 key_name: " + (key + "." + key_child).lower())
                                    if compare_name.lower() == (key + "." + key_child).replace("_", "").lower():
                                        statName = ".".join(name_list) + "." + key + "." + key_child
                                        if isinstance(value_child, dict):
                                            statName = statName + "." + list(value_child.items())[0][0]
                                            value_child = list(value_child.items())[0][1]
                                        if "countMin" in statName:
                                            statName = statName.replace("countMin", "").replace("_", "")
                                        complete_dict[name_count].append(value_child)
                                        complete_dict[name_count].append("")
                                        complete_dict[name_count].append("")
                                        # print("        " + statName + " => " + value_child)
                        else:
                            name = name.replace('cpu0', 'cpu')
                            if name not in complete_dict:
                                complete_dict[name] = []
                                complete_dict[name].append(name)
                            complete_dict[name].append(v)
                            complete_dict[name].append(value)
                            complete_dict[name].append("")
                            complete_dict[name].append("")
                            # print(name + " => " + v)
                            # print("        " + ".".join(name_list) + "." + key + " => " + value)


def compare(path):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Comparing stats file in {path}')  # Press ⌘F8 to toggle the breakpoint.
    files = os.listdir(path)
    if "stats.txt" not in files:
        print(f"stats.txt not found in {path}")
        return
    stats = path + "/stats.txt"
    statsRead = open(stats, "r")
    dataList = {}
    for line in statsRead:
        if "#" not in line:
            continue
        dataSplit = line.split("#")
        stat = list(filter(lambda empty: empty != "", dataSplit[0].split(" ")))
        statName = stat[0].strip()
        statCount = stat[1].strip()
        statDesc = dataSplit[1].strip()
        # print(f"{statDesc}: {statName} => {statCount}")
        statNameComponents = statName.split(".")
        parsedList = stat_name_to_list(statNameComponents, statCount, 0)
        dict_merge(dataList, parsedList)

    for k, v in dataList.items():
        find_similar_count_min(k, v, v, [k])
    statsRead.close()
    return


def compare_all(directory):
    global directory_
    strategy = ['nc', 'c', 'm', 'mc']
    directorye = directory.split("/")[0] + "/"
    for benchmark in list_benchmarks:
        bench = benchmark.split(".")[1]
        dpath = directorye + bench + "/"
        direc = dpath + "compared/"
        if 'compared' not in os.listdir(dpath):
            os.mkdir(direc)
        for s in size:
            for d in depth:
                for i in range(0, 5):
                    for strat in strategy:
                        directory_ = dpath + str(s) + "_" + str(d) + "_" + strat + "_" + str(i)
                        try:
                            compare(directory_)
                        except FileNotFoundError:
                            print(directory_ + " not found")

                    comparePathFile = direc + str(s) + "_" + str(d) + "_" + str(i) + "_compared.csv"
                    compareFile = open(comparePathFile, "w")
                    csvWriter = csv.writer(compareFile, "excel")
                    print(f"-----------------------------------------------------\n"
                          f"Writing to file: {comparePathFile}\n"
                          f"-----------------------------------------------------")
                    for key in complete_dict:
                        csvWrite = False
                        for index in range(1, (len(complete_dict[key]) - 1), 2):
                            if complete_dict[key][index] != '0' and complete_dict[key][index] != '':
                                csvWrite = True
                        if csvWrite:
                            # print(complete_dict[key])
                            for x in range(0, 17 - len(complete_dict[key])):
                                complete_dict[key].append("")
                            csvWriter.writerow(complete_dict[key])
                    complete_dict.clear()


filter_files_individual = []


def filter_all(f):
    global filter_files_individual
    filter_files_individual = []
    counter_list = [1, 2, 3, 4, 0]

    directorye = f
    if f[-1] != "/":
        directorye = directorye + ""
    filterFile = open("filter.txt", "r")
    filter_list = []
    for la in filterFile:
        filter_list.append(la.replace("\n", ""))
    direc = directorye
    files = os.listdir(direc)
    if f"filtered" not in files:
        print(f"Create Directory {direc}/filtered")
        os.mkdir(f"{direc}/filtered")
    for s in size:
        for d in depth:
            all_data = []
            for counter in counter_list:
                file = f"{s}_{d}_{counter}_compared.csv"
                if file not in files:
                    print(f"{file} not found in {directory}")
                    continue

                csvFile = f"{direc}{file}"
                print(csvFile)
                compareFile = open(csvFile, "r")
                csvReader = csv.reader(compareFile, "excel")
                for line in csvReader:
                    if line[0] in filter_list:
                        all_data.append(line)
                all_data.append(['\n'])
                compareFile.close()
            if len(all_data):
                filteredFile = f"{direc}filtered/{s}_{d}_filtered.csv"
                filterFile = open(filteredFile, "w")
                filter_files_individual.append(filteredFile)
                csvWriter = csv.writer(filterFile, "excel")
                print(f"Writing to file {filteredFile}")
                for line in all_data:
                    csvWriter.writerow(line)
                filterFile.close()


def combine_filter_files(f):
    global filter_files_individual
    filesData = []
    directorye = f

    if f[-1] != "/":
        directorye = directorye + "/"

    for fileName in filter_files_individual:
        openF = open(fileName, "r")
        csvReader = csv.reader(openF)
        data = []
        i = 0
        data.append([])
        data[i] = {}
        for line in csvReader:
            if len(line) > 1:
                # a = 1
                data[i][line[0].replace('cpu0', 'cpu')] = line
            else:
                i += 1
                data.append([])
                data[i] = {}
        filesData.append(data)

    finalData = [{}, {}, {}, {}, {}]

    posIndex = 0
    for group in filesData:
        groupIndex = 0
        for mode in group:
            for key, event in mode.items():
                if key not in finalData[groupIndex].keys():
                    finalData[groupIndex][key] = []
                # if len(finalData[posIndex][key]) != 0 and len(finalData[posIndex][key]) < groupIndex * 17:
                if len(finalData[groupIndex][key]) < posIndex * 17:
                    # empty = [""] * 17
                    # finalData[groupIndex][key] += empty
                    finalData[groupIndex][key] += event
                else:
                    finalData[groupIndex][key] += event
            groupIndex += 1
        posIndex += 1

    if len(finalData):
        filteredFile = f"{directorye}filtered/filtered.csv"
        filterFile = open(filteredFile, "w")
        filter_files_individual.append(filteredFile)
        csvWriter = csv.writer(filterFile, "excel", delimiter=';')
        print(f"Writing to file {filteredFile}")
        modes = ['cms', 'cms-cu', 'cms-morris', 'morris']

        dataTable = [
            ['32_2', '256_2'],
            ['64_2', '512_2'],
            ['64_4', '512_4'],
            ['128_4', '1024_4'],
            ['128_8', '1024_8'],
            ['256_4', '2048_4'],
            ['256_8', '2048_8'],
            ['512_4', '4096_4'],
            ['512_8', '4096_8'],
            ['512_16', '4096_16'],
            ['1024_4', '8192_4'],
            ['1024_8', '8192_8'],
            ['1024_16', '8192_16'],
        ]

        for group in finalData:
            csvWriter.writerow(excel_labels * 12)
            min_val = []
            max_val = []
            rmsd_val = []
            nrmsd_val = []
            len_events = []
            combinedData = {}
            for pos in range(0, 13):
                combinedData[pos] = {}
                for i in range(0, 4):
                    combinedData[pos][modes[i]] = {
                        'real_values': [],
                        'sim_values': [],
                        'sums': 0,
                    }
                    for key in group:
                        if len(group[key]) >= (pos * 17 + i * 4 + 1) and group[key][pos * 17 + i * 4 + 1] != '':
                            # if key == 'system.cpu.lsq0.countMinSquashedLoads':
                            # print(i == 3 and int(group[key][pos*17+i*4+1]) == 0)
                            if i == 3 and int(group[key][pos * 17 + i * 4 + 1]) == 0:
                                continue
                            combinedData[pos][modes[i]]['real_values'].append(int(group[key][pos * 17 + i * 4 + 2]))
                            combinedData[pos][modes[i]]['sim_values'].append(int(group[key][pos * 17 + i * 4 + 1]))
                            group[key][pos * 17 + i * 4 + 3] = int(group[key][pos * 17 + i * 4 + 2]) - int(
                                group[key][pos * 17 + i * 4 + 1])
                            group[key][pos * 17 + i * 4 + 4] = group[key][pos * 17 + i * 4 + 3] * group[key][
                                pos * 17 + i * 4 + 3]
                            combinedData[pos][modes[i]]['sums'] += group[key][pos * 17 + i * 4 + 4]
            for key, lines in group.items():
                # if key == 'system.cpu.lsq0.countMinSquashedLoads':
                # print(lines)
                csvWriter.writerow(lines)
            for pos in range(0, 13):
                min_val.append('Min')
                max_val.append('Max')
                rmsd_val.append('RMSD')
                nrmsd_val.append('NRMSD')
                len_events.append('NUM_EVENTS')
                for i in range(0, 4):
                    min_val.append('')
                    max_val.append('')
                    rmsd_val.append('')
                    nrmsd_val.append('')
                    len_events.append(len(combinedData[pos][modes[i]]['real_values']))
                    dataTable[pos].append(str(len(combinedData[pos][modes[i]]['real_values'])))

                    minimum = ''
                    if len(combinedData[pos][modes[i]]['real_values']) > 0:
                        minimum = (min(combinedData[pos][modes[i]]['real_values']))
                    min_val.append(minimum)
                    maximum = ''
                    if len(combinedData[pos][modes[i]]['real_values']) > 0:
                        maximum = (max(combinedData[pos][modes[i]]['real_values']))
                    max_val.append(maximum)
                    rmsd_val.append('')
                    nrmsd_val.append('')
                    len_events.append('')

                    min_val.append('')
                    diff = ''
                    if maximum != '' and minimum != '':
                        diff = maximum - minimum
                    max_val.append(str(diff))
                    rmsd_val.append('')
                    nrmsd_val.append('')
                    len_events.append('')

                    sumValues = combinedData[pos][modes[i]]['sums']
                    averageValue = ''
                    rmsd = ''
                    if len(combinedData[pos][modes[i]]['real_values']) > 0:
                        averageValue = sumValues / len(combinedData[pos][modes[i]]['real_values'])
                        rmsd = float(math.sqrt(averageValue))
                    min_val.append(averageValue)
                    max_val.append('')

                    rmsd_val.append(str(rmsd).replace('.', ','))
                    nrmsd = ''
                    if rmsd != '':
                        nrmsd = float(rmsd / (maximum - minimum))
                    nrmsd_val.append(str(nrmsd).replace('.', ','))
                    len_events.append('')
                    dataTable[pos].append(str(nrmsd).replace('.', ','))
                dataTable[pos] += ['******', dataTable[pos][0], dataTable[pos][1]]
            csvWriter.writerow([""])
            csvWriter.writerow([""])
            csvWriter.writerow([""])
            csvWriter.writerow([""])
            csvWriter.writerow([""])
            csvWriter.writerow(min_val)
            csvWriter.writerow(max_val)
            csvWriter.writerow(rmsd_val)
            csvWriter.writerow(nrmsd_val)
            csvWriter.writerow(len_events)
        for dataT in dataTable:
            csvWriter.writerow(dataT)
        filterFile.close()


plot_save = True


def plot_seperate(f, img_bench_folder, plotType):
    filtered = f + "filtered.csv"
    if 'filtered.csv' not in os.listdir(f):
        print(f"Filtered.csv not found in {filtered}")
        return
    csvFile = open(filtered, "r")
    csvReader = csv.reader(csvFile, delimiter=";")
    csvData = []
    for i in csvReader:
        csvData.append(i)
    data = {}
    ignore = [-11, -9, -7, -5, -4, -2, -1]
    for i in range(-13, 0):
        if i in ignore:
            continue
        line = csvData[i]
        gr = 1
        for index in range(0, 5):
            size_index = line[index * 11 + 1]  # .split("_")[0]
            if size_index not in data.keys():
                data[size_index] = {}
            if str(gr) not in data[size_index].keys():
                data[size_index][str(gr)] = {}
            size_group = data[size_index][str(gr)]
            if 'cms' not in data[size_index][str(gr)]:
                size_group['cms'] = [[int(line[index * 11 + 2].replace(",", "."))],
                                     [float(line[index * 11 + 3].replace(",", "."))]]
            else:
                size_group['cms'][0].append(int(line[index * 11 + 2].replace(",", ".")))
                size_group['cms'][1].append(float(line[index * 11 + 3].replace(",", ".")))
            if 'cms-cu' not in size_group:
                size_group['cms-cu'] = [[int(line[index * 11 + 4].replace(",", "."))],
                                        [float(line[index * 11 + 5].replace(",", "."))]]
            else:
                size_group['cms-cu'][0].append(int(line[index * 11 + 4].replace(",", ".")))
                size_group['cms-cu'][1].append(float(line[index * 11 + 5].replace(",", ".")))
            if 'cms-morris' not in size_group:
                size_group['cms-morris'] = [[int(line[index * 11 + 6].replace(",", "."))],
                                            [float(line[index * 11 + 7].replace(",", "."))]]
            else:
                size_group['cms-morris'][0].append(int(line[index * 11 + 6].replace(",", ".")))
                size_group['cms-morris'][1].append(float(line[index * 11 + 7].replace(",", ".")))
            if 'morris' not in size_group:
                size_group['morris'] = [[int(line[index * 11 + 8].replace(",", "."))],
                                        [float(line[index * 11 + 9].replace(",", "."))]]
            else:
                size_group['morris'][0].append(int(line[index * 11 + 8].replace(",", ".")))
                size_group['morris'][1].append(float(line[index * 11 + 9].replace(",", ".")))
            gr = gr + 1

    if plotType == 'events':
        saveEventsNRMSD(data, img_bench_folder)
    else:
        saveSizeNRMSD(data, img_bench_folder)
    return


allDataAverage = {}


def saveSizeNRMSD(data, img_bench_folder):
    global allDataAverage
    averaged = {}
    bench = img_bench_folder.split("/")[-1]
    allDataAverage[bench] = {}
    for si, values in data.items():
        siz = si.split("_")[0]
        for gr, strat_vals in values.items():
            if gr not in averaged.keys():
                averaged[gr] = {}
            if siz not in averaged[gr].keys():
                averaged[gr][siz] = {}

            if gr not in allDataAverage[bench].keys():
                allDataAverage[bench][gr] = {}
            if siz not in allDataAverage[bench][gr].keys():
                allDataAverage[bench][gr][siz] = {}
            for strategy, valu in strat_vals.items():
                if strategy not in averaged[gr][siz].keys():
                    averaged[gr][siz][strategy] = [[], []]
                if strategy not in allDataAverage[bench][gr][siz].keys():
                    allDataAverage[bench][gr][siz][strategy] = [[], []]

                averaged[gr][siz][strategy][0] = averaged[gr][siz][strategy][0] + [
                    int(round(sum(data[si][gr][strategy][0]) / len(data[si][gr][strategy][0]), 0))]
                averaged[gr][siz][strategy][1] = averaged[gr][siz][strategy][1] + [
                    sum(data[si][gr][strategy][1]) / len(data[si][gr][strategy][1])]

                allDataAverage[bench][gr][siz][strategy][0] = allDataAverage[bench][gr][siz][strategy][0] + [
                    int(round(sum(data[si][gr][strategy][0]) / len(data[si][gr][strategy][0]), 0))]
                allDataAverage[bench][gr][siz][strategy][1] = allDataAverage[bench][gr][siz][strategy][1] + [
                    float(sum(data[si][gr][strategy][1]) / len(data[si][gr][strategy][1]))]

    avg = {}
    for gr, values in averaged.items():
        if gr not in avg.keys():
            avg[gr] = {}
        for si, strat_vals in values.items():
            if si not in avg[gr].keys():
                avg[gr][si] = {}
            for strategy, valu in strat_vals.items():
                avg[gr][si][strategy] = [
                    int(round(sum(averaged[gr][si][strategy][0]) / len(averaged[gr][si][strategy][0]), 0)),
                    sum(averaged[gr][si][strategy][1]) / len(averaged[gr][si][strategy][1])]
            # print(avg[gr][si])
    plotSizeNRMSD(avg, img_bench_folder)


def saveEventsNRMSD(data, img_bench_folder):
    global allDataAverage
    averaged = {}
    bench = img_bench_folder.split("/")[-1]
    allDataAverage[bench] = {}

    for si, values in data.items():
        siz = si.split("_")[0]
        if siz not in averaged.keys():
            averaged[siz] = {}
        if siz not in allDataAverage[bench].keys():
            allDataAverage[bench][siz] = {}
        for gr, strat_vals in values.items():
            if gr not in averaged[siz].keys():
                averaged[siz][gr] = {}
            if gr not in allDataAverage[bench][siz].keys():
                allDataAverage[bench][siz][gr] = {}
            for strategy, valu in strat_vals.items():
                if strategy not in averaged[siz][gr].keys():
                    averaged[siz][gr][strategy] = [[], []]
                if strategy not in allDataAverage[bench][siz][gr].keys():
                    allDataAverage[bench][siz][gr][strategy] = [[], []]

                averaged[siz][gr][strategy][0] = averaged[siz][gr][strategy][0] + [
                    int(round(sum(data[si][gr][strategy][0]) / len(data[si][gr][strategy][0]), 0))]
                averaged[siz][gr][strategy][1] = averaged[siz][gr][strategy][1] + [
                    float(sum(data[si][gr][strategy][1]) / len(data[si][gr][strategy][1]))]

                allDataAverage[bench][siz][gr][strategy][0] = allDataAverage[bench][siz][gr][strategy][0] + [
                    int(round(sum(data[si][gr][strategy][0]) / len(data[si][gr][strategy][0]), 0))]
                allDataAverage[bench][siz][gr][strategy][1] = allDataAverage[bench][siz][gr][strategy][1] + [
                    float(sum(data[si][gr][strategy][1]) / len(data[si][gr][strategy][1]))]

    avg = {}
    for si, values in averaged.items():
        if si not in avg.keys():
            avg[si] = {}
        for gr, strat_vals in values.items():
            if gr not in avg[si].keys():
                avg[si][gr] = {}
            for strategy, valu in strat_vals.items():
                avg[si][gr][strategy] = [
                    int(round(sum(averaged[si][gr][strategy][0]) / len(averaged[si][gr][strategy][0]), 0)),
                    sum(averaged[si][gr][strategy][1]) / len(averaged[si][gr][strategy][1])]
    plotEventsNRMSD(avg, img_bench_folder)


def allSizeNRMSD(img_bench_folder):
    averaged = {}
    for bench, dataAverage in allDataAverage.items():
        for gr, values in dataAverage.items():
            for si, strat_vals in values.items():
                if gr not in averaged.keys():
                    averaged[gr] = {}
                if si not in averaged[gr].keys():
                    averaged[gr][si] = {}

                for strategy, valu in strat_vals.items():
                    if strategy not in averaged[gr][si].keys():
                        averaged[gr][si][strategy] = [[], []]
                    averaged[gr][si][strategy] = [
                        int(round(sum(allDataAverage[bench][gr][si][strategy][0]) / len(
                            allDataAverage[bench][gr][si][strategy][0]), 0)),
                        sum(allDataAverage[bench][gr][si][strategy][1]) / len(
                            allDataAverage[bench][gr][si][strategy][1])]
    plotSizeNRMSD(averaged, img_bench_folder)


def allEventsNRMSD(img_bench_folder):
    averaged = {}
    for bench, dataAverage in allDataAverage.items():
        for si, values in dataAverage.items():
            if si not in averaged.keys():
                averaged[si] = {}
            for gr, strat_vals in values.items():
                if gr not in averaged[si].keys():
                    averaged[si][gr] = {}

                for strategy, valu in strat_vals.items():
                    if strategy not in averaged[si][gr].keys():
                        averaged[si][gr][strategy] = [[], []]
                    averaged[si][gr][strategy] = [
                        int(round(sum(allDataAverage[bench][si][gr][strategy][0]) / len(
                            allDataAverage[bench][si][gr][strategy][0]), 0)),
                        sum(allDataAverage[bench][si][gr][strategy][1]) / len(
                            allDataAverage[bench][si][gr][strategy][1])]
    plotEventsNRMSD(averaged, img_bench_folder)


def plotSizeNRMSD(averaged, img_bench_folder):
    plots = {}
    strategy = ['cms', 'cms-cu', 'cms-morris', 'morris']
    marker = {
        'cms': "v",
        'cms-cu': "^",
        'cms-morris': "s",
        'morris': "x"
    }

    for gr, gr_strategy in averaged.items():
        # print(gr)
        for s, strateg in gr_strategy.items():
            if gr not in plots.keys():
                plots[gr] = {}
            for strat, values in strateg.items():
                if strat not in plots[gr].keys():
                    plots[gr][strat] = {}
                    plots[gr][strat]['size'] = [int(s)]
                    plots[gr][strat]['NRMSD'] = [float(averaged[gr][s][strat][1])]
                    plots[gr][strat]['num_events'] = [averaged[gr][s][strat][0]]
                else:
                    plots[gr][strat]['size'].append(int(s))
                    plots[gr][strat]['NRMSD'].append(float(averaged[gr][s][strat][1]))
                    plots[gr][strat]['num_events'].append(averaged[gr][s][strat][0])

        plt.clf()
        num_events = []
        maxNRMSD = 0
        minNRMSD = 1000
        for strat in strategy:
            if strat == 'morris':
                continue
            m = max(plots[gr][strat]['NRMSD'])
            mi = max(plots[gr][strat]['NRMSD'])
            if m > maxNRMSD:
                maxNRMSD = m
            if mi < minNRMSD:
                minNRMSD = mi
            cmsx = np.array(plots[gr][strat]['size'])
            cmsy = np.array(plots[gr][strat]['NRMSD'])
            # print(plots[gr][strat]['size'])
            # print(plots[gr][strat]['NRMSD'])
            plt.plot(cmsx, cmsy, marker=marker[strat], markersize=7)
            if strat != 'morris':
                num_events = (num_events + plots[gr][strat]['num_events'])
            else:
                i = 0
                add = 0.1
                for s in plots[gr][strat]['size']:
                    # print((int(s), plots[gr][strat]['NRMSD'][i]))
                    # print((int(s), float(plots[gr][strat]['NRMSD'][i]) + add * (i + 1.4)))
                    # print()
                    y = float(plots[gr][strat]['NRMSD'][i]) + add * (i - 1.4)
                    if y < 0.03:
                        y = 0.03
                    plt.annotate(str(plots[gr][strat]['num_events'][i]),
                                 xy=(int(s), plots[gr][strat]['NRMSD'][i]),
                                 xytext=(int(s), y),
                                 arrowprops=dict(facecolor='black', shrink=1, width=0.1, headwidth=5, headlength=5)
                                 )
                    i = i + 1
        limy = maxNRMSD * 1.2
        if limy > 1:
            limy = 1
        num_events = int(round(sum(num_events) / len(num_events), 0))
        plt.plot([(num_events * 64), (num_events * 64)], [0, maxNRMSD], marker='.')

        plt.suptitle(f"Size VS NRMSD for ~{num_events} events")
        plt.xlabel("Size in Bits")
        plt.ylabel("NRMSD")
        plt.legend(strategy[0:-1] + ['baseline'])
        plt.xscale('log', base=2)
        # plt.yscale('log')
        plt.gca().xaxis.set_major_formatter(tck.ScalarFormatter())
        plt.gca().yaxis.set_major_formatter(tck.ScalarFormatter())
        plt.gca().yaxis.set_major_locator(tck.LinearLocator(7))
        plt.gca().set_xlim(450, 12000)
        plt.gca().set_ylim(0, limy)
        formatImg = "png"
        img = f"{img_bench_folder}/size_nrmsd_{num_events}.{formatImg}"
        if f"{num_events}.{formatImg}" in os.listdir(img_bench_folder):
            os.remove(img)

        plt.grid(True, which='both')
        if plot_save:
            plt.savefig(img, dpi=512, format=formatImg)
        else:
            plt.show()
    return


def plotEventsNRMSD(averaged, img_bench_folder):
    plots = {}
    strategy = ['cms', 'cms-cu', 'cms-morris', 'morris']
    marker = {
        'cms': "v",
        'cms-cu': "^",
        'cms-morris': "s",
        'morris': "x"
    }

    for si, gr_strategy in averaged.items():
        # print(gr)
        for gr, strateg in gr_strategy.items():
            if si not in plots.keys():
                plots[si] = {}
            for strat, values in strateg.items():
                if strat not in plots[si].keys():
                    plots[si][strat] = {}
                    plots[si][strat]['gr'] = [int(gr)]
                    plots[si][strat]['NRMSD'] = [float(averaged[si][gr][strat][1])]
                    plots[si][strat]['num_events'] = [averaged[si][gr][strat][0]]
                else:
                    plots[si][strat]['gr'].append(int(gr))
                    plots[si][strat]['NRMSD'].append(float(averaged[si][gr][strat][1]))
                    plots[si][strat]['num_events'].append(averaged[si][gr][strat][0])

        plt.clf()
        siz = []
        maxNRMSD = 0
        minNRMSD = 1000
        for strat in strategy:
            m = max(plots[si][strat]['NRMSD'])
            mi = max(plots[si][strat]['NRMSD'])
            if m > maxNRMSD:
                maxNRMSD = m
            if mi < minNRMSD:
                minNRMSD = mi
            cmsx = np.array(plots[si][strat]['num_events'])
            cmsy = np.array(plots[si][strat]['NRMSD'])
            # print(plots[si][strat]['num_events'])
            # print(plots[si][strat]['NRMSD'])
            plt.plot(cmsx, cmsy, marker=marker[strat], markersize=7)
            # if strat != 'morris':
            #    si = si
            siz = (siz + [int(si)])
            # else:
            #     i = 0
            #     add = 0.1
            #     for s in plots[si][strat]['size']:
            #         # print((int(s), plots[gr][strat]['NRMSD'][i]))
            #         # print((int(s), float(plots[gr][strat]['NRMSD'][i]) + add * (i + 1.4)))
            #         # print()
            #         y = float(plots[si][strat]['NRMSD'][i]) + add * (i - 1.4)
            #         if y < 0.03:
            #             y = 0.03
            #         plt.annotate(str(plots[si][strat]['num_events'][i]),
            #                      xy=(int(s), plots[si][strat]['NRMSD'][i]),
            #                      xytext=(int(s), y),
            #                      arrowprops=dict(facecolor='black', shrink=1, width=0.1, headwidth=5, headlength=5)
            #                      )
            #         i = i + 1

        limy = maxNRMSD * 1.2
        if limy > 1:
            limy = 1
        size_str = int(sum(siz) / len(siz))
        # num_events = int(round(sum(num_events) / len(num_events), 0))
        plt.plot([(size_str / 64), (size_str / 64)], [0, limy * 0.7], marker='.')
        #
        # print(siz)
        plt.suptitle(f"Number of Events VS NRMSD for {size_str} Bits")
        plt.xlabel("Number of Events")
        plt.ylabel("NRMSD")
        plt.legend(strategy + ['baseline'])
        plt.xscale('log', base=2)
        # plt.yscale('log')
        plt.gca().xaxis.set_major_formatter(tck.ScalarFormatter())
        plt.gca().yaxis.set_major_formatter(tck.ScalarFormatter())
        plt.gca().yaxis.set_major_locator(tck.LinearLocator(7))
        plt.gca().set_xlim(10, 192)
        plt.gca().set_ylim(0, limy)
        formatImg = "png"
        img = f"{img_bench_folder}/events_nrmsd_{size_str}.{formatImg}"
        if f"{size_str}.{formatImg}" in os.listdir(img_bench_folder):
            os.remove(img)

        plt.grid(True, which='both')
        if plot_save:
            plt.savefig(img, dpi=512, format=formatImg)
        else:
            plt.show()
    return

def printHelp():
    print("Usage --- python3 main.py {command} {directory} {plotType for plot command}")
    print("command = {'compare', 'filter', 'plot'}")
    print("directory = Directory filled with gem5 stats")
    print("plotType = {'events', 'size'}")

if __name__ == '__main__':
    print("Analysing Stats")
    if len(sys.argv) <= 2:
        if len(sys.argv) <= 1:
            print("No command provided")
            printHelp()
        else:
            print("No directory/File provided for comparing stats")
            printHelp()
        exit(1)
    operation = sys.argv[1]
    directory = sys.argv[2]
    if operation == "compare":
        print("Comparing Stats")
        compare_all(directory)
    elif operation == "filter":
        print("Filtering Stats")
        for bm in list_benchmarks:
            b = bm.split(".")[1]
            folder = directory.split("/")[0] + "/" + b + "/compared/"
            filter_all(folder)
            combine_filter_files(folder)
    elif operation == "plot":
        print("Plotting Stats")
        plotType = ""
        if len(sys.argv) > 3:
            plotType = sys.argv[3]
        if len(sys.argv) > 4:
            plot_save = sys.argv[4]
            if plot_save == "show":
                plot_save = False
        if plotType in ['events', 'size']:
            print("undefined PlotType")
            exit(1)
        if "img" not in os.listdir("./"):
            os.mkdir("img")
        if "individual" not in os.listdir("./img"):
            os.mkdir("img/individual")
        for bm in list_benchmarks:
            b = bm.split(".")[1]
            print(f"Plotting graph for {b}")
            if b not in os.listdir("./img/individual"):
                os.mkdir(f"./img/individual/{b}")
            folder = directory.split("/")[0] + "/" + b + "/compared/filtered/"
            plot_seperate(folder, f"./img/individual/{b}", plotType)

        if "mean" not in os.listdir("./img"):
            os.mkdir("./img/mean")
        if plotType == 'events':
            allEventsNRMSD("./img/mean")
        else:
            allSizeNRMSD("./img/mean")
    else:
        printHelp()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
