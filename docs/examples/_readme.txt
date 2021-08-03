To create a new example, follow these steps:

1. Create a new python script for the example in /examples/. You should be
   able to run this script by itself, and produce a plot at the end.
2. Create a .rst file in this directory that has:
   - A title
   - .. jupyter-execute:: setup.py
        :hide-code:

   - .. jupyter-execute:: ../../examples/FILENAME.py


   where FILENAME is the name of the file created in step 1.
