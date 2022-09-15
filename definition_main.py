import json
import glob
import os
import csv
import re

#Load .json file
with open("input_data.json","r") as f:
    json_data = json.load(f)

#Get every value from .json file
front_elements = json_data["front_elements"]
back_elements = json_data["back_elements"]
spaces = json_data["spaces"]
alignment_name = json_data["alignment_name"]
input_file_path = json_data["input_file_path"]
output_file_path = json_data["output_file_path"]

def confirm_files():
    #Confirm input file
    print("-----------------")
    txt_files=glob.glob(r'{}/*.txt'.format(input_file_path)) #LIST
    fasta_files=glob.glob(r'{}/*.fasta'.format(input_file_path)) #LIST
    search_data =txt_files+fasta_files #LIST
    print(search_data)
    print("-----------------")
    return search_data

def info_print(*args):
    print("-----------------")
    for i in args:
        print(str(i))
    print("-----------------")

def logic_and_output(name,sequence,N_count):
    #Initialize
    cordence_rows = []
    semi_cordence_rows = []
    header = ["front_elements","spaces",
        "back_elements","start_num","end_num","region"
        ] #columns: 1~6
    os.makedirs(output_file_path, exist_ok=True)

    #The Length of back_elements, front_elements
    b_len = len(back_elements)
    f_len = len(front_elements)

    #Find "front_elements" from "sequence"
    iterator_front = re.finditer(front_elements,sequence) #ITERATOR
    if any(iterator_front): #True/False
        for i in iterator_front:
            group = i.group()
            s = i.start()+N_count
            e = i.end()+N_count #end_num
            region =str(s)+":"+str(e)
            print("\t"+group+str(s)+str(e))
            try:
                expected_b = sequence[e+spaces:e+spaces+b_len]
                if expected_b == back_elements: #cordence
                    print("\t"+"  [cordence sequence EXISTS!]")
                    cordence_rows.append([front_elements,spaces,back_elements,s,e,region])
                else: #semi_cordence
                    print("\t"+"  Only semi_cordence")
                    semi_cordence_rows.append([front_elements,spaces,expected_b,s,e,region])
            except IndexError as IE:
                print(IE)
            except PermissionError as PE:
                print(PE)
    else:
        print("\t"+str(front_elements)+" is not found")

    #Find "back_elements" from "sequence"
    iterator_back = re.finditer(back_elements,sequence) #ITERATOR
    if any(iterator_back):#True/False
        for i in iterator_back:
            group = i.group()
            s = i.start()+N_count #start_num
            e = i.end()+N_count
            region =str(s)+":"+str(e)
            print("\t"+group+str(s)+str(e))
            try:
                expected_f = sequence[s-spaces-f_len:s-spaces]
                if expected_f == front_elements:#cordence
                    print("\t"+"  [cordence sequence EXISTS!]")
                    cordence_rows.append([front_elements,spaces,back_elements,s,e,region]) #columns: 1~6
                else:#semi_cordence
                    print("\t"+"  Only semi_cordence")
                    semi_cordence_rows.append([expected_f,spaces,back_elements,s,e,region]) #columns: 1~6
            except IndexError as IE:
                print(IE)
            except PermissionError as PE:
                print(PE)
    else:
        print("\t"+str(back_elements)+" is not found")

    #Write "cordence_rows" on .csv file
    cordence_="./cordence_result_"+name+".csv"
    cordence_csv=r"{}".format(output_file_path+cordence_)
    with open(cordence_csv,"w",newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(cordence_rows)

    #Write "semi_cordence_rows" on .csv file
    semi_cordence_="./semi_cordence_result_"+name+".csv"
    semi_cordence_csv=r"{}".format(output_file_path+semi_cordence_)
    with open(semi_cordence_csv,"w",newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(semi_cordence_rows)


def main():
    search_data = confirm_files() #LIST
    total = str(len(search_data))
    print("=================START=================")
    for idx,search_file in enumerate(search_data,start=1):
        with open(search_file,"r") as f1:
            #Initialize
            name = f1.readline()
            num = 5 if search_file.endswith(".txt") else 1
            name = name[num:40] #info_print(arg), logic_and_output(arg)

            #Initialize AND Pass args to info_print()
            sequence=f1.read() #logic_and_output(arg)
            N_count=sequence.count("N")
            N_info = "N_counts: "+str(N_count) #info_print(arg)
            base_pairs = str(len(sequence))+"bp" #info_print(arg)
            sequence=sequence.replace("N", "").replace(" ", "").replace("\n", "") #info_print(arg)
            info_print(name,N_info,base_pairs)

            #Search
            logic_and_output(name,sequence,N_count)
        print("===============continue("+str(idx)+"/"+total+")==============") #Example: ==continue(8/24)==
    print("[All completed]: "+total+" files!")
    print("==================END==================")

if __name__ == "__main__":
    main()
