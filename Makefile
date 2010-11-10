ROOT=/data/cygob2

install: imgview.py
	cp imgview.py $(ROOT)/share/imgview/

install_db: datastore.db3
	cp datastore.db3 $(ROOT)/tables/
