# InsertFaviconLink.py

from os import *
from os.path import *
from re import *

print "Copying the favicon to the lib directory...";
system("cp img/favicon.ico lib/");

# insert the link to the favicon before </html>
print "Inserting the link to the favicon...";

f = file("lib/index.html", 'r')          # load file into a string
text = f.read()
f.close()

text = text.replace("<title>",'<link rel="shortcut icon" href="favicon.ico" type="image/vnd.microsoft.icon" />\n<title>')      
                                                        
f = file("lib/index.html", 'w')          # write file again
f.write(text)
f.close()    

print "ready!\n"
