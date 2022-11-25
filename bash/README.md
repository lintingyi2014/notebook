
### Bash scripting  
UNIX/ Linux  

Declaring a variable looks like this:  
`variable="Some string"` (no spaces)

Using the variable:  
```
echo "$variable" # => Some string
echo '$variable' # => $variable
```
If you want to use the variable's value, you should use $.

Parameter expansion ${...}:  
`echo "${variable}" # => Some string`  
This is a simple usage of parameter expansion such as two examples above.
Parameter expansion gets a value from a variable.
It "expands" or prints the value.
During the expansion time the value or parameter can be modified.
Below are other modifications that add onto this expansion.
