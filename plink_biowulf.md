
### Using plink1/2 interchangeably in Biowulf


The python script says: 
```
    parser.add_argument('--plink1_path', type=str, default='plink',
                        help='Path to plink. Needed if plink is not '
                             'in $PATH')
    parser.add_argument('--plink2_path', type=str, default='plink2',
                        help='Path to plink. Needed if plink is not '
                             'in $PATH')
```

So one way around it is to put the 2 plink directories into your path. In your .bashrc, add
```
export PATH=/usr/local/apps/plink/3.6-alpha:/usr/local/apps/plink/1.9.0-beta4.4:$PATH
```
#re-read the updated bashrc file
```
source ~/.bashrc
```
### check that plink and plink2 are in the path
```
echo $PATH  
```
### you should see the directories in your path
```
which plink
```
### the command above should report
/usr/local/apps/plink/3.6-alpha/plink2
```
which plink2
```
### the command above should report
/usr/local/apps/plink/3.6-alpha/plink2

