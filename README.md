Grade: A3;  To implement a chemical compound database in Python and SQL, make a graphical interface that allows the user to interact and filter compounds in the database. 

Input Required:
1. A collection of 100 compounds in .sdf file format was used as input data.

What the code does:
1. It uses various modules from rdkit to extract and/or calculate 14 fields important for evaluating the compound’s suitability to be a hit.
2. A relational database was created using the SQLite3 module in Python.
3. To handle this database, a class called DatabaseManager was created, that allows the creation of the database with the ‘CREATE’ query
 provided, loads the database using the ‘INSERT’ query provided and finally returns the search of a ‘SELECT’ statement provided.
4.  The user interface was programmed using the Tkinter module in Python. The database from SQLite3 is queried every time any data has to be loaded into the interface.
5.  Many features are provided to the user including, searching for a compound name, using preset rules like Lipinski's Rule of Five to filter the compounds or
 even manually specify values to filter the compounds.
6. Finally, the filtered compounds can be saved as a .csv file for further analysis.
