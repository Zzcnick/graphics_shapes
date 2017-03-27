all: Picture.java
	javac Picture.java && java Picture

run: Picture.java
	make && make png && echo "Saved as out.png" && display out.png

clean: 
	rm *.class *.ppm *~ *.jpg *.png

jpg: out.ppm
	java Picture; \
	convert out.ppm out.jpg

png: out.ppm
	java Picture; \
	convert out.ppm out.png
