.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<

phy=gbritoap@phycalc2.physique.ens-lyon.fr

folder=/home/gbritoap/SPDE-JKW-Algorithm/

## send : send .py and .sh files from local to phycalc2
.PHONY: send
send:
	rsync -avR *.py $(phy):$(folder)
	rsync -avR *.sh $(phy):$(folder)
	rsync -avR Makefile $(phy):$(folder)

## recv : recv results from phycalc2 to local
.PHONY: recv
recv:
	rsync -avr $(phy):$(folder)/results/* results/
	rsync -avr $(phy):$(folder)/data/*_R_*000_N_* data/
