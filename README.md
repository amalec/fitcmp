**fitcmp.py** is written by Adrian Malec. This program is free software (MIT
licensed).

# LICENSE

Copyright (C) 2012 Adrian Malec

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# NOTE

I've used some bits and pieces of code from the Internet. I've indicated where
and have included the URL. I use Pexpect Copyright (c) 2008 Noah Spurrier.

# ABOUT

This python script displays vpfit-like plot for two fort.13 files for 
comparison. Since I first wrote it at the start of my PhD I've added further
features to the script. For more information use the -h flag or just run the
script. The code is not my best Python code, but is the best I could do given
very limited time I could spend on this.

The dependencies are as follows:
* Python 2.7.2
* matplotlib 1.1.0
* numpy 1.6.1
* asciitable 0.8

You may be able to run the script with lower versions of these. Maybe not.

# RUNNING FITCMP.PY

You run the script by executing fitcmp (the main script). You may like to add
its location to your $PATH.

You must always run the script with an argument specifying the path of
your vpfit executable. Creating an alias to your vpfit in your working 
directory may be handy.

Examples:

* fitcmp ./vpfit -o old_fort.13 -n fort.13
* fitcmp ./vpfit -n fort.13/initial4.fit

Your vpfit must be correctly configured, including environmental paths such as
ATOMDIR, VPFSETUP and VPFPLOTS. This script uses the ATOMDIR in your
environment.

Note that you may easily customise many aspects of the output plots by 
creating a matplotlibrc file in your current working directory.
See: http://matplotlib.sourceforge.net/users/customizing.html
This may be useful for quickly creating figures for publications or 
presentations. In this case, experiment with figure.dpi as your starting
point.

# DOCUMENTATION

Formal documentation is nonexistent, sorry. Most operations are
self-explanatory and are listed, along with brief descriptions, by typing '?'
at the input prompt.

# SUPPORT

Please report bugs and problems at https://github.com/amalec/fitcmp/
I personally provide no guarantees or support outside those granted by the 
goodness of my heart.

I would appreciate your acknowledgement if you find this script useful, but
it is of course not necessary. I would also be grateful if you sent me any 
major changes or bug fixes, so I can implement them in my copy.
