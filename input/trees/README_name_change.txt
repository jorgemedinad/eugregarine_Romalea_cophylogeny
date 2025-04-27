# Name change.

In Windows.

Open terminal.

Make sure you have your name change file (name_change.txt) is in the same folder as your tree.

You need a tree file (yourtree.tre) which has only the specimen codes.

# The format of the name_change.txt is as follows. You should have the list matching the species list.
s/specimen_code/whatever_you_want_to_change_to/

#In terminal, type the following. Open Ubuntu in Windows. Go to working directory using format cd /mnt/c/Users/path/to/file
# Do not tip names and do not use parenthesis
sed -f name_change.txt yourtree.tre > yourtree_renamed.tre