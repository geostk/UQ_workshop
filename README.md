# UQ_workshop
Scripts and related materials for Uncertainty in Geo-science. A Workshop on Hazard Analysis Buffalo March 15-18 2016

## Instruction to use
1. Download the repository as a zip file from here:  <a href="https://github.com/haghakhani/UQ_workshop/archive/master.zip" target=0>UQ_workshop</a> 
2. Log in to your account at VHUB and open a new workspace
3. Upload the zip file in your space at VHUB
4. Unzip the file by running: 
    <code>unzip UQ_workshop-master.zip -d .</code> 
5. Go to the directory by running:
    <code>cd UQ_workshop-master/</code> 
6. The next step is to install numpy library. First extract the source code folder by running: 
        <code>tar -xvf numpy-1.10.4.tar.gz</code>
7. Then got to the extracted folder by typing the following in the terminal:
    <code>cd numpy-1.10.4/</code>
8. To install numpy you need to run the following command:
    <code>python setup.py install --user</code> 
9. Now go back to the main directory:
     <code>cd ..</code> 
10. Now you need to add the new installed library into your python path by:
     <code>setenv PYTHONPATH $PWD:$PYTHONPATH</code> 
11. Now you should be able to run the scripts as described in their instructions inside the python or ipython
    
      
