MSG="Changes"
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

branch:
	@-echo $(BRANCH)

clean:
	@-find . -name "*~" -exec rm -rf {} \;
	@-find . -name "#*#" -exec rm -rf {} \;
	@-find . -name "*.pyc" -exec rm -rf {} \;
	@-find . -name ".DS_Store" -exec rm -rf {} \;

cleandata:
	@-find . -name "__pycache__" -exec rm -rf {} \; >& /dev/null
	@-rm -rf confs/* tmp/* plots/*

cleanall:clean cleandata

commit:cleanall
	@-git commit -am "$(MSG)"
	@-git push origin $(BRANCH)

pull:
	git reset --hard HEAD	
	git pull
