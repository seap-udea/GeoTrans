clean:
	-find . -name "*~" -exec rm -rf {} \;
	-find . -name "#*#" -exec rm -rf {} \;
	-find . -name "*.pyc" -exec rm -rf {} \;

cleandata:
	-find . -name "__pycache__" -exec rm -rf {} \;
	rm -rf confs/* tmp/* plots/*

cleanall:clean cleandata

commit:
	-git commit -am "Commit"
	git push origin master

pull:
	git reset --hard HEAD	
	git pull
