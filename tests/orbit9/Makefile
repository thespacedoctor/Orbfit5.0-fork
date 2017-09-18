
tests: testbaric testprop testo9

testbaric:
	mv inbaric.opt inbaric.opt.tmp
	cp -p testout/inbaric.opt .
	time inbaric.x
	time inbarmerc.x
	-(diff planxx.inc testout/)
	-(diff pvenexx.inc testout/)
	-(diff barsunxx.inc testout/)
	-(diff barmerxx.inc testout/)
	mv inbaric.opt.tmp inbaric.opt

testprop:
	cp -p testout/prop.opt .
	time prop9.x
	-(diff v9.mea testout/v9.mea)
	-(diff v9.rem testout/v9.rem)
	-(diff v9.pro testout/v9.pro)
	-(diff v9.ang testout/v9.ang)
	-(diff v9.dis testout/v9.dis)
	-(diff v9.adi testout/v9.adi)
	-(diff v9.res testout/v9.res)

testo9:
	mv orb9.opt orb9.opt.tmp
	cp -p testout/orb9.opt .
	time orbit9.x
	-(diff vast.dat testout/vast.dat)
	-(diff vpla.dat testout/vpla.dat)
	-(diff vast.fil testout/vast.fil)
	-(diff vpla.fil testout/vpla.fil)
	-(diff orb9.out testout/orb9.out)
	mv orb9.opt.tmp orb9.opt

testo9t:
	time inbaric.x
	time orbit9t.x
	echo 'v' | time propert3.x

clean:
	-(rm orb9.dma orb9.dmp orb9.num orb9.out v*.fil v*.dat orb9.clo)
	-(rm v9* giffv.gnu giffv.tmp out.* gmon.out *~ *.done)
	-(rm barmerxx.inc pvenexx.inc allplmxx.inc barsunxx.inc allplaxx.inc planxx.inc epsi.out)
	-(rm vprt* libcenter.dat orb9t.sum orb9t.out)

distclean: clean