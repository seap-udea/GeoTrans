MSG="Changes"
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

branch:
	@-echo $(BRANCH)

clean:
	@-find . -name "*~" -exec rm -rf {} \;
	@-find . -name "#*#" -exec rm -rf {} \;
	@-find . -name "*.pyc" -exec rm -rf {} \;

cleandata:
	@-find . -name "__pycache__" -exec rm -rf {} \;
	@-rm -rf confs/* tmp/* plots/*

cleanall:clean cleandata

commit:
	@-git commit -am "$(MSG)"
	@-git push origin $(BRANCH)

pull:
	git reset --hard HEAD	
	git pull
