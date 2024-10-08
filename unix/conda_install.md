# Conda Linux installation instructions

Conda can be installed for a single user, so each of you will install your own copy on the server.
You can do this on any computer, including your own, which is why we want you to know about it.
"Mini"conda simply means the kit mostly just contains conda, without any extra stuff.

First, log on to the course server.

From your browser, get the link to the installation kit from the official site:
https://docs.conda.io/en/latest/miniconda.html#linux-installers

Right-click > copy link on the "Miniconda3 Linux 64-bit".

In the terminal, write `wget` and paste the link after it. The press enter.
The command should look like:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Now make the kit executable. In Unix systems, the file extension is actually quite irrelevant.
It's the file properties and content that determine how they're handled by the system.
The `chmod` command changes file properties. 
Below, `u` means user (that's you!) and `+x` means "set executable" (for the user only).

```
chmod u+x Miniconda3-latest-Linux-x86_64.sh
```

To run an executable from the current directory, run `./EXECUTABLE_NAME` (no space after `./`).
The dot means "this directory" since Unix ignores the current directory by default when looking for executables.

```
./Miniconda3-latest-Linux-x86_64.sh
```



The installation is interactive:

- It will ask you to read and accept a license agreement.
Just pres [Enter] until it asks you to accept by writing "yes".

- It will then tell you where conda will be installed. 
The default location is fine, so just press [Enter]

- The next question is "Do you wish the installer to initialize Miniconda3
by running conda init? [yes|no]". Write "yes" (no quotes) and press [Enter].
The visible effect of this is that when you'll open a new terminal, you'll see
`(base)` to the left of the command prompt. This just means conda is "active" for you.

if `(base)` is not there, run the following command
```
source ~./bashrc
```


That's it!

To save space, remove the kit :)




```
rm Miniconda3-latest-Linux-x86_64.sh
```

# Troubleshooting 

If for some reason conda is not initalized (i.e. you don't see the `(base)` marker),
run `source ~/.bashrc` in the terminal.



