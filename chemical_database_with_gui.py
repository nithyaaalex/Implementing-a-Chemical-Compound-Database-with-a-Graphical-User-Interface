#import statements
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
import sqlite3
import os
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
import tkinter.simpledialog
from ttkbootstrap import Style
import csv


#declaring some global variables
path_working_directory = "D:/coding_assignment/"
molecules = []
photo_images = []
db_file = "2933044_Chemical_Database.db"
createdb_file = "sql_create_commands.txt"
input_database_file = "Molecules17.sdf"
columns_filter = {'MolecularWeight': '4', 'logD': '5', 'logP': '6', 'HBondAcceptors': '7', 'HBondDonors': '8', 'RingCount': '9', 'AromaticFusedRings': '10', 'PolarSurfaceArea': '11', 'RotatableBonds': '12'}
filter_entries = {}
#setting the working directory
os.chdir(path_working_directory)


#class definition for database handling
class Database_Manager:
    def __init__(self, db_file):
        """
        Initialises the object with the filename of the sqlite database 
        """
        self.database_filename = db_file
    
    def createdb(self, query):
        """
        This function executes one or multiple CREATE statements using executescript on the database specified
        Keyword arguments:
        query -- can be one or multiple CREATE statements
        """
        connection = sqlite3.connect(self.database_filename)
        cursor = connection.cursor()
        error  = False
        #if no errors occur during creation of tables/database, the changes will be committed, else they'll be rolled back
        try:
            cursor.executescript(query)
        except Exception as e:
            connection.rollback()
            print(f"An exception has been encountered {e} and the create statements (DDL) did not complete")
            error = True
        if error == False:
            connection.commit()
        cursor.close()
        connection.close()

    def loaddb(self, query, parameters):
        """
        This function executes a group of INSERT queries on the database specified
        Keyword arguments:
        query -- INSERT statement 
        parameters -- list of (tuples of parameters) for the INSERT statement
        """
        connection = sqlite3.connect(self.database_filename)
        cursor = connection.cursor()
        error = False
        #if no errors occur during insertion, the changes will be committed, else they'll be rolled back
        try:
            cursor.executemany(query, parameters) #need to supply list of parameters as tuples
        except Exception as e:
            connection.rollback()
            print(f"An exception has been encountered {e} and the insert query with parameters did not complete")
            error = True
        if error == False:
            connection.commit()
        cursor.close()
        connection.close()
    
    def querydb(self, query):
        """
        This function executes a single sql SELECT query on the database specified and returns the result of the transaction
        It returns an iterable object with each line as a row retrived from the database
        Keyword arguments:
        query -- select statement 
        """
        connection = sqlite3.connect(self.database_filename)
        cursor = connection.cursor()
        rows = cursor.execute(query).fetchall()
        connection.commit()
        cursor.close()
        connection.close()
        return rows


#MAIN starts here
#parsing and calculating fields from the sdf file & storing them into a list
with Chem.MultithreadedSDMolSupplier(input_database_file) as sdf_file_contents:
  for mol in sdf_file_contents: #going through each row in the sdf file
    if mol is not None:
        #for each row, extracting data from the sdf (name, id, formula and logd are extracted directly)
        id = "Compound_" + mol.GetProp("CdId")
        name = mol.GetProp("Name")
        formula = mol.GetProp("Formula")
        log_D = mol.GetProp("LogD")

        #creates the images if not already present
        if not os.path.exists(f"image{id}.png"):
            img = Draw.MolToImage(mol)
            img.save(f"image{id}.png")

        #all the following parameters were calculated using rdkit 
        smiles_string = Chem.MolToSmiles(mol)
        molecular_weight = '{:.2f}'.format(Descriptors.MolWt(mol))
        log_P = '{:.2f}'.format(Descriptors.MolLogP(mol))
        h_bond_acceptors = Descriptors.NumHAcceptors(mol)
        h_bond_donors = Descriptors.NumHDonors(mol)
        ring_count = Descriptors.RingCount(mol)  
        aromatic_fused_rings = sum(1 for ring in Chem.GetSymmSSSR(mol) if len(ring) > 1 and all(mol.GetBondWithIdx(bond_idx).GetIsAromatic() for bond_idx in ring))        
        polar_surface_area = '{:.2f}'.format(Descriptors.TPSA(mol))
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)

        #saving all the values as a tuple in a list
        molecules.append((id, name, formula, smiles_string, molecular_weight, log_D, log_P, h_bond_acceptors, h_bond_donors, ring_count, aromatic_fused_rings, polar_surface_area, rotatable_bonds))

#creating an object of the database manager class to interact with it
omics_database = Database_Manager(db_file) 
#creating and loading the database if it doesn't exist in the file path
if not os.path.exists(db_file):
    with open(createdb_file) as create_table_query:
        create_query = create_table_query.read() #storing all the lines as one big string
    omics_database.createdb(create_query) 
    general_query = 'INSERT INTO ChemicalDatabase VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)' 
    omics_database.loaddb(general_query, molecules)


#function definitions 
    
#this function deletes all the rows in the view and loads new rows in according to the query supplied
def loading_data(query_input):
    query_output = omics_database.querydb(query_input)
    chemical_treeview.delete(*chemical_treeview.get_children())
    i= 0
    for row in query_output:
        #formatting the image to input it into table
        ID_for_i = row[0]
        tk_image = Image.open(f"image{ID_for_i}.png")
        tk_image.thumbnail((100,100))
        tk_image = ImageTk.PhotoImage(tk_image)
        i+=1
        photo_images.append(tk_image)
        #inserting the image as well as all other rows together
        chemical_treeview.insert(parent = '', index = 'end', image=tk_image, values=row)
    count_label.config(text=f"Total Rows Showing: {i}") #changes the count 
    #adds the rows back in
    chemical_treeview.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

#default for loading the whole database back in
def load_all():
    query_input = "SELECT * FROM ChemicalDatabase;"
    loading_data(query_input=query_input)

#each of the following functions has a different query according to the filter needed, and calls loading data to show the result
def lipinski_filter(): 
    #displays a pop-up dialog 
    tkinter.simpledialog.messagebox.showinfo("Parameters used:", "MolecularWeight <= 500, logP <= 5, HBondAcceptors <= 10, HBondDonors <= 5")
    query_input = "SELECT * FROM ChemicalDatabase WHERE (MolecularWeight <= 500 AND logP <= 5 AND HBondAcceptors <= 10 AND HBondDonors <= 5);"
    loading_data(query_input=query_input)

def lead_likeness_filter(): 
    #displays a pop-up dialog 
    tkinter.simpledialog.messagebox.showinfo("Parameters used:", "MolecularWeight <= 450, (logD <= 4 AND logD >= -4), RingCount <= 4, RotatableBonds <= 10, HBondAcceptors <= 8, HBondDonors <= 5")
    query_input = "SELECT * FROM ChemicalDatabase WHERE (MolecularWeight <= 450 AND (logD <= 4 AND logD >= -4) AND RingCount <= 4 AND RotatableBonds <= 10 AND HBondAcceptors <= 8 AND HBondDonors <= 5);"
    loading_data(query_input=query_input)

def bioavailability_filter(): 
    #displays a pop-up dialog
    tkinter.simpledialog.messagebox.showinfo("Parameters used:", "any 6 conditions true out of : MolecularWeight <= 500, logP <= 5, HBondAcceptors <= 10, HBondDonors <= 5, RotatableBonds <= 10, PolarSurfaceArea <= 200, FusedAromaticRings <= 5")
    query_input = "SELECT * FROM ChemicalDatabase WHERE (CASE WHEN MolecularWeight <= 500 THEN 1 ELSE 0 END + CASE WHEN logP <= 5 THEN 1 ELSE 0 END + CASE WHEN HBondDonors <= 5 THEN 1 ELSE 0 END + CASE WHEN HBondAcceptors <= 10 THEN 1 ELSE 0 END + CASE WHEN RotatableBonds <= 10 THEN 1 ELSE 0 END + CASE WHEN PolarSurfaceArea <= 200 THEN 1 ELSE 0 END + CASE WHEN FusedAromaticRings <= 5 THEN 1 ELSE 0 END )>= 6";
    loading_data(query_input=query_input)

#function for the search bar - gets the user input and searches the database for all matches (even partial)
def search():
    search_term = search_entry.get()
    loading_data(query_input="SELECT * FROM ChemicalDatabase WHERE CommonName LIKE '%" + search_term + "%'")

#filters data based on individual columns (min and max given by user)
def filter_data():
    #base sql query, will build on this 
    query = "SELECT * FROM ChemicalDatabase WHERE "
    conditions = [] #only those conditions that have atleast one non zero value will be filtered
    for col, index in columns_filter.items():
        min_val = float(filter_entries[index]['min_entry'].get())
        max_val = float(filter_entries[index]['max_entry'].get())
        #according to which value has been entered by the user, a differeny query will be appended to the main
        #if this is not done, it will create queries like greater than 100 but less than 0
        if min_val == 0 and max_val != 0:
            conditions.append(f"{col} <= {max_val}") 
        elif min_val != 0 and max_val == 0:
            conditions.append(f"{col} >= {min_val}")
        elif min_val != 0 and max_val != 0:
            conditions.append(f"{col} BETWEEN {min_val} AND {max_val}")
    #if there is atleast one condition supplied, then only it will proceed
    if conditions:
        if len(conditions) == 1:
            query += conditions[0]  #only one condition, no need for AND in sql
        else:
            query += ' AND '.join(conditions) #adds all the conditions with AND inbetween
        print(query) #for verifying the query if needed
        loading_data(query_input=query)
    
#to order by based on the dropdown
def sort_table_based_on(value):
    query = "SELECT * FROM ChemicalDatabase ORDER BY " + value
    loading_data(query_input=query)

#saves the current rows of the ui into a csv, allows user to choose name of the file
def save_csv_with_dialog():
    #asks user for filename
    file_name = tkinter.simpledialog.askstring("Save As", "Enter file name:")
    if file_name:
        #adds the csv extension if not provided by the user
        if not file_name.endswith(".csv"):
            file_name += ".csv"
        with open(file_name, "w", newline='') as myfile:
            csvwriter = csv.writer(myfile, delimiter=',')
            csvwriter.writerow(columns)
            for row_id in chemical_treeview.get_children():
                row = chemical_treeview.item(row_id)['values']
                csvwriter.writerow(row)
        #tells the user once it has been saved
        tkinter.messagebox.showinfo("Save Successful", f"File '{file_name}' saved successfully.")


#GUI part
#creates the main window
main_window = tk.Tk()
main_window.title('Chemical Databases - 2933044')

#using bootstrap preset styles
style = Style(theme='vapor')
style.configure("Treeview",rowheight=104)

#creates a frame for title
title_frame = ttk.Frame(main_window)
title_frame.pack(side=tk.TOP, padx=10, pady=10)

#frame for filters
filter_frame = ttk.Frame(main_window)
filter_frame.pack(side=tk.RIGHT, padx=10, pady=10)

#frame for table and scrollbar
table_frame = ttk.Frame(main_window)
table_frame.pack(fill=tk.BOTH, padx=10, pady=10)

#creates vertical scrollbar
scroll_y = ttk.Scrollbar(table_frame, orient='vertical')
scroll_y.pack(side=tk.RIGHT, fill=tk.Y)

#creates horizontal scrollbar
scroll_x = ttk.Scrollbar(table_frame, orient='horizontal')
scroll_x.pack(side=tk.BOTTOM, fill=tk.X)

#making columns and adding scrollbar
columns = ('ID', 'Name', 'Formula', 'SMILES', 'Molecular Weight', 'Log D','Log P', 'H Bond Acceptors', 'H Bond Donors', 'Ring Count', 'Aromatic Fused Rings', 'Polar Surface Area', 'Rotatable Bonds')
chemical_treeview = ttk.Treeview(table_frame, columns=columns,yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set) 

#configures the scrollbar
scroll_y.config(command=chemical_treeview.yview)
scroll_x.config(command=chemical_treeview.xview)

#centering and changing all column widths
for col in columns:
    chemical_treeview.heading(col, text=col, anchor=tk.CENTER)
    chemical_treeview.column(col, anchor=tk.CENTER, width=150)  

#to center and increase width for image column, and label that column
chemical_treeview.column("#0", width=150)
chemical_treeview.heading("#0", anchor=tk.CENTER, text='Compound Image')

#title for whole thing
title_label = ttk.Label(title_frame, text="Chemical Structure Databases", font=("Helvetica", 18))
title_label.grid(row=0, column=2, padx=5, pady=5) 

#label that will display the count of rows
count_label = ttk.Label(title_frame, text="Total Rows Showing: 0",font=("Helvetica", 14))
count_label.grid(row=2, column=2, padx=5, pady=5)

#search bar stuff, labels and entry widget
search_label = ttk.Label(title_frame, text="Enter Compound Name To Search:")
search_label.grid(row=3, column=2, padx=5, pady=5) 
#creates an entry widget for the search query
search_entry = ttk.Entry(title_frame, width=30)
search_entry.grid(row=4, column=2, padx=5, pady=5, sticky="ew", columnspan=2)

#binds the enter key to the search function (so user can click on search or hit enter and it will work)
search_entry.bind('<Return>', lambda event: search())

#creates a button to trigger the search
search_button = ttk.Button(title_frame, text="Search", command=search,bootstyle="danger")
search_button.grid(row=4, column=4, padx=5, pady=5, sticky="e")

#loading all 100 entries in
load_all()

#creating some frame for some labels
label_frame_top= ttk.Frame(filter_frame)
label_frame_top.grid(row=0, column=1, padx=5, pady=5, sticky='n')
#creating some more labels
label_rules = ttk.Label(label_frame_top, text="Preset Rule Filter:")
label_rules.pack()

#all the filter buttons
filter_button1 = ttk.Button(filter_frame, text="Lipinski's Rule", command = lipinski_filter,bootstyle="secondary")
filter_button1.grid(row=1, column=0, padx=20, pady=20, sticky='n')
filter_button2 = ttk.Button(filter_frame, text="Lead-Likeness", command = lead_likeness_filter,bootstyle="warning")
filter_button2.grid(row=1, column=1, padx=20, pady=20, sticky='n')
filter_button3 = ttk.Button(filter_frame, text="Bioavailability", command = bioavailability_filter,bootstyle="info")
filter_button3.grid(row=1, column=2, padx=20, pady=20, sticky='n')

#ORDER by label 
order_entry_label = ttk.Label(filter_frame, text="Order By:")
order_entry_label.grid(row=2, column=0, padx=5, pady=5, sticky='n')    

#DROPDOWN for orderby 
#dropdown_options = ["MolecularWeight","logD","logP","HBondAcceptors","HBondDonors","RingCount","FusedAromaticRings", "PolarSurfaceArea","RotatableBonds"]
dropdown_options = ["MolecularWeight","MolecularWeight DESC", "logD","logD DESC","logP","logP DESC", "HBondAcceptors","HBondAcceptors DESC", "HBondDonors","HBondDonors DESC" ,"RingCount","RingCount DESC", "FusedAromaticRings","FusedAromaticRings DESC", "PolarSurfaceArea","PolarSurfaceArea DESC","RotatableBonds", "RotatableBonds DESC"]
#to store the options for dropdown
selected_option = tk.StringVar()

#setting the initial value that is shown off
selected_option.set(dropdown_options[0])

#actual menu option widget
dropdown_menu = tk.OptionMenu(filter_frame, selected_option, *dropdown_options, command=sort_table_based_on)
dropdown_menu.grid(row=2, column=1, padx=5, pady=5, sticky='n')

#another frame for some labels
label_frame2 = ttk.Frame(filter_frame)
label_frame2.grid(row=4, column=0, padx=5, pady=5, sticky='n')
#creating some labels
label_filter = ttk.Label(label_frame2, text="Individual Row Filter:")
label_filter.pack()

#labels to show which entrybox is which (for min and max of the filter part)
min_label = ttk.Label(filter_frame, text="Min") 
min_label.grid(row=len(filter_entries)+4, column=1, sticky='n')  
max_label = ttk.Label(filter_frame, text="Max")
max_label.grid(row=len(filter_entries)+4, column=2, sticky='n')  

#creates filter entry fields for each column
i=5 #for positionining purposes
for col, index in columns_filter.items():
    #writing the column name (the name of the value to be filtered on)
    label = ttk.Label(filter_frame, text=f"{col}:")
    label.grid(row=i, column=0, padx=5, pady=5, sticky='n')

    #creates an entry bar for the minimum value
    min_entry = ttk.Entry(filter_frame)
    min_entry.grid(row=i, column=1, padx=5, pady=5, sticky='n')
    min_entry.insert(tk.END, '0.0')  #SETTING default to zero

    #creates entry for max value
    max_entry = ttk.Entry(filter_frame)
    max_entry.grid(row=i, column=2, padx=5, pady=5, sticky='n')
    max_entry.insert(tk.END, '0.0')  #SETTING default to zero
    i+=1

    #saves the values in a list, as a dictionary, by index
    filter_entries[index] = {'min_entry': min_entry, 'max_entry': max_entry}

#creates a button to apply the filter (actually executes)
filter_button = ttk.Button(filter_frame, text="Apply Filter", command=filter_data)
filter_button.grid(row=i, column=1, padx=20, pady=20, sticky='n')

#more buttons - save and show all buttons (extra features)
save_button =ttk.Button(filter_frame,text="Save",command=save_csv_with_dialog)
save_button.grid(row=i, column=2, padx=20, pady=20, sticky='n')
filter_button0= ttk.Button(filter_frame, text="Show All", command = load_all)
filter_button0.grid(row=i, column=0, padx=20, pady=20, sticky='n')


#another frame for to position a label for general instruction to user
label_frame = ttk.Frame(filter_frame)
label_frame.grid(row=i+1, column=1, padx=5, pady=5, sticky='n')
#creating some labels
label = ttk.Label(label_frame, text="Click on 'Show All' to see all compounds again")
label.pack()

#starts the main event loop
main_window.mainloop()

#COULD HAVE IMPLEMENTED
#add argument parser to allow user to change database name, path, and the other global variables using flags
#could have implemented the sort on the records that are visible in the gui not on the whole database
#logger to log each event in the UI (clicks, saves, filters, searches)
#search bar dynamic (as user is typing, the compounds that match will show up as a suggestion (like google autocomplete))
