GitHub for Noobs
================

Make a GitHub repository from a local folder
--------------------------------------------

1. On your GitHub profile, click the topleft `+ > New repository`. Name it as you please (*e.g.* your local folder's name for clarity).

**Do not Initialise with a README**, we want to initialise with the local folder.

*(Conversely, if you wish to start your project directly on GitHub, then initialising with a README is definitely the thing to do.)*

Click on `Clone with HTTPS` and simply copy the web URL to clipboard.

*****

2. In a console, `cd` to the local folder.

Run

	git init
*Intilialise the local folder as a git repository, creating a .git folder where changes will be tracked.*

	git remote add origin <paste_https_link>
*The local folder is now linked to the GitHub created in 1.*

	git add .
*Add the content of the folder to git.*

	git commit -m "Initial commit"
*Log the change. The `-m` option allows you to add a message to keep track of things clearly.*

	git push origin master
*Put everything new/modified online on GitHub.*

Git may ask for your identity at some point, in which case just comply running

	git config user.email "you@example.com"
	git config user.name "Your Name"

Now, you should be set up both locally on the original machine and online on GitHub.


Using GitHub: good practices
----------------------------

To download for the first time the repository on a new machine, do

    git clone <https_link_to_github_repository>

Do this only the first time, as it will download **everything**.

Start every day with

	cd path/to/local/folder
	git pull origin master

It will update your local files with any update made anywhere else.

Whenever you edit/add a file, wherever you do it (on your professional machine, on your personal one, online on GitHub, ...), when you're done don't forget to

	git commit -m "Tell about what you did."
	git push origin master

Did you forget whether you were productive today? `git status` will remind you!
