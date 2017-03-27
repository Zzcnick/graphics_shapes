import java.io.*;
import java.util.*;

public class Picture {

    public static void main(String[] args) throws FileNotFoundException {
	if (args.length > 0) { // Parser
	    Canvas c = new Canvas(500, 500, 255, 255, 255);
	    String script = args[0];
	    File file = new File(script);	    
	    Scanner sc = new Scanner(file);

	    String cmd = "";
	    while (sc.hasNext()) {
		cmd = sc.next(); // Command
		
		// Processing
		if (cmd.equals("line")) {
		    c.edge(sc.nextDouble(), sc.nextDouble(), sc.nextDouble(),
			   sc.nextDouble(), sc.nextDouble(), sc.nextDouble());
		} else if (cmd.equals("bezier")) {
		    c.bezier(sc.nextDouble(), sc.nextDouble(), sc.nextDouble(), sc.nextDouble(),
			     sc.nextDouble(), sc.nextDouble(), sc.nextDouble(), sc.nextDouble());
		} else if (cmd.equals("hermite")) {
		    c.hermite(sc.nextDouble(), sc.nextDouble(), sc.nextDouble(), sc.nextDouble(),
			      sc.nextDouble(), sc.nextDouble(), sc.nextDouble(), sc.nextDouble());
		} else if (cmd.equals("circle")) {
		    c.circle(sc.nextDouble(), sc.nextDouble(), sc.nextDouble(), sc.nextDouble());
		} else if (cmd.equals("scale")) {
		    c.scale(sc.nextDouble(), sc.nextDouble(), sc.nextDouble());
		} else if (cmd.equals("move")) {
		    c.translate(sc.nextDouble(), sc.nextDouble(), sc.nextDouble());
		} else if (cmd.equals("rotate")) {
		    c.rotate(sc.next().charAt(0), sc.nextDouble());
		} else if (cmd.equals("apply")) {
		    c.apply();
		} else if (cmd.equals("ident")) {
		    c.clearTransform();
		} else if (cmd.equals("draw")) {
		    c.draw();
		} else if (cmd.equals("save")) {
		    c.save(sc.next());
		}
	    }
	    return;
	} 
	// Curves  ============================
	Canvas c = new Canvas(500, 500, 0, 40, 60);
	Pixel p = new Pixel(240, 250, 255);
	// int frame = 0;

 	// c.savestate(); // Background
	for (int i = 0; i < 3; i++) {
	    // Left Eye
	    c.hermite(100,325+i,200,325+i,50+i,75+i,50+i,-75+i,p);
	    c.circle(160, 325, 0, 18+i, p);
	    
	    // Right Eye
	    c.hermite(300,325+i,400,325+i,50+i,75+i,50+i,-75+i,p);
	    c.circle(360, 325, 0, 18+i, p);
	    
	    // Nose
	    c.line(230 + i, 340,
		   230 + i, 290, p);
	    c.bezier(230 + i, 290 - i, 305 + i, 320 + i,
		     305 + i, 210 - i, 230 + i, 240 + i, p);
		
	    // Mouth
	    c.hermite(170,220-i,330,220-i,
		      75, -100, 75, 100, p);
	}
	c.draw();
	c.save("out.ppm");

	// GIF
	// for (int r = 0; r < 360; r += 3) {
	//     c.load();
	//     frame++;

	//     String buffer = "" + frame;
	//     while (buffer.length() < 3)
	//    	buffer = "0" + buffer;
	//     c.draw();
	//     c.save(buffer + ".ppm");
	    
	//     c.translate(-250,0,0);
	//     c.rotate('y', 3);
	//     c.translate(250,0,0);
	//     c.apply();
	// }
	// ==================================== */
	
	/* // Transform ==========================
	   Canvas c = new Canvas(500, 500, 20, 0, 20);
	   double x = 0; double y = 0;
	   double r = 0; double a = 0;
	   int R, G, B; R = 0; G = 0; B = 0;
	   double tmp = 25 * Math.random();
	   int frame = 0;

	   // Background
	   for (int i = 0; i < 150; i++) {
	   R = 155 + (int)(100 * Math.random());
	   G = 155 + (int)(100 * Math.random());
	   B = 155 + (int)(100 * Math.random());
	   Pixel p = new Pixel(R, G, B);
	   r = 2 + (int)(3 * Math.random());
	   x = 20 + (int)(460 * Math.random());
	   y = 20 + (int)(460 * Math.random());
	   c.edge(x - r, y, x + r, y, p);
	   c.edge(x, y - r, x, y + r, p);
	   }
	   c.draw();
	   c.clearEdges();
	   c.savestate();

	   // Inner Ring - Back
	   for (int i = 0; i < 250; i++) {
	   R = 180 + (int)(50 * Math.random());
	   G = R;
	   B = 20;
	   Pixel p = new Pixel(R, G, B);
	   r = 50 + 50 * Math.random();
	   a = Math.random() * Math.PI - Math.PI / 2;
	   x = r * Math.cos(a);
	   y = r * Math.sin(a);
	   double extension = 5 + 10 * Math.random();
	   double dx, dy;
	   dx = extension * Math.sin(a);
	   dy = extension * Math.cos(a);
	   c.edge(x + dx,
	   y - dy,
	   x - dx,
	   y + dy,
	   p);
	   }

	   // Nucleus
	   for (int fold = 0; fold < 360; fold += 10) {
	   for (int i = 0; i < 360; i += (int)(25 * Math.random())) {
	   R = 150 + (int)(106 * Math.random());
	   G = 100;
	   B = 20;
	   Pixel p = new Pixel(R, G, B);
	   x = 40 * Math.cos(2 * Math.PI * i / 360);
	   y = 40 * Math.sin(2 * Math.PI * i / 360);
	   c.edge(0,0,x,y,p);
	   }
	   c.rotate('x', 10);
	   c.apply();
	   }

	   // Inner Ring - Front
	   for (int i = 0; i < 250; i++) {
	   R = 180 + (int)(50 * Math.random());
	   G = R;
	   B = 30;
	   Pixel p = new Pixel(R, G, B);
	   r = 50 + 50 * Math.random();
	   a = Math.random() * Math.PI + Math.PI / 2;
	   x = r * Math.cos(a);
	   y = r * Math.sin(a);
	   double extension = 5 + 10 * Math.random();
	   double dx, dy;
	   dx = extension * Math.sin(a);
	   dy = extension * Math.cos(a);
	   c.edge(x + dx,
	   y - dy,
	   x - dx,
	   y + dy,
	   p);
	   }

	   // Outer Ring
	   for (int i = 0; i < 250; i++) {
	   R = 130 + (int)(50 * Math.random());
	   G = R - 10;
	   B = 10;
	   Pixel p = new Pixel(R, G, B);
	   r = 100 + 100 * Math.random();
	   a = Math.random() * 2 * Math.PI;
	   x = r * Math.cos(a);
	   y = r * Math.sin(a);
	   double extension = 15 + 10 * Math.random();
	   double dx, dy;
	   dx = extension * Math.sin(a);
	   dy = extension * Math.cos(a);
	   c.edge(x + dx,
	   y - dy,
	   x - dx,
	   y + dy,
	   p);
	   }

	   c.rotate('x', 60);
	   c.rotate('y', 30);
	
	   // while (frame < 180) {
	   //    String buffer = "" + frame;
	   //    while (buffer.length() < 3)
	   // 	buffer = "0" + buffer;

	   // c.rotate('y', 2);
	   c.translate(250, 250, 250);
	   c.apply();

	   c.draw();

	   c.save("out.ppm");

	   // c.translate(-250, -250, -250);
	   // c.apply();

	   // c.load();
	   // frame++;
	   // } // End Animation Loop

	   // ==================================== */

	/* // EdgeMatrix =========================
	   Canvas c = new Canvas(500, 500, 0, 0, 0);
	   double x, y, X, Y;
	   Pixel p = new Pixel();
	   boolean odd = false;
	   x = 250; y = 250;

	   for (int i = 0; i < 720; i++) {
	   if (i % 90 == 0) { // Move Pivot
	   x = 250 + 100 * Math.cos(2 * Math.PI * i / 720.);
	   y = 250 + 100 * Math.sin(2 * Math.PI * i / 720.);
	   i += 10;
	   if (i % 4 == 0) {
	   p = new Pixel(255, 128, 10); // Orange
	   odd = true;
	   } else {
	   p = new Pixel(10, 128, 255); // Blue
	   odd = false;
	   }
	   }
	   X = 250 + 225 * Math.cos(2 * Math.PI * i / 720. + Math.PI / 7.8);
	   Y = 250 + 225 * Math.sin(2 * Math.PI * i / 720. + Math.PI / 7.8);
	   if (odd)
	   c.edge(x, y, X, Y, new Pixel(p.adjust(-2, -1, 0)));
	   else
	   c.edge(x, y, X, Y, new Pixel(p.adjust(0, -1, -2)));
	   }

	   c.draw();
	   c.save("out.ppm");
	   // ==================================== */

	/* // Lines ==============================
	   Canvas c = new Canvas(500, 500, 0, 0, 0);

	   Pixel p = new Pixel(18, 10, 143);
	   int x1, y1, x2, y2;
	   for (int i = 0; i < 500; i += 5) {
	   x1 = i; y1 = 499;
	   x2 = 499; y2 = 499 - i;
	   c.line(x1, y1, x2, y2, p);
	   }
	   for (int i = 0; i < 500; i += 10) {
	   x1 = 499; y1 = 499 - i;
	   x2 = 499 - i; y2 = 0;
	   c.line(x1, y1, x2, y2, p);
	   }
	   for (int i = 0; i < 500; i += 15) {
	   x1 = 499 - i; y1 = 0;
	   x2 = 0; y2 = i;
	   c.line(x1, y1, x2, y2, p);
	   }
	   for (int i = 0; i < 500; i += 20) {
	   x1 = 0; y1 = i;
	   x2 = i; y2 = 499;
	   c.line(x1, y1, x2, y2, p);
	   }
	   c.save("out.ppm");
	   // ==================================== */

	/* // Line Testing  ===================
	   c = new Canvas();
	   c.line(250, 250, 400, 250); // +x
	   c.line(250, 250, 350, 300); // 1
	   c.line(250, 250, 300, 350); // 2
	   c.line(250, 250, 250, 400); // +y
	   c.line(250, 250, 200, 350); // 3
	   c.line(250, 250, 150, 300); // 4
	   c.line(250, 250, 100, 250); // -x
	   c.line(250, 250, 150, 200); // 5
	   c.line(250, 250, 200, 150); // 6
	   c.line(250, 250, 250, 100); // -y
	   c.line(250, 250, 300, 150); // 7
	   c.line(250, 250, 350, 200); // 8
	   // ==================================== */
	
	/* // Green Mountains =================
	   Canvas c = new Canvas(1000, 500, 16, 0, 32);
	   for (int y = 480; y > 0; y -= 10) {
	   for (int x = 0; x < 1000; x++) {
	   if (Math.random() < 0.005) { // 0.5%
	   int g = (int)(256 * (1 - y / 500.));
	   int r = (int)(g * 0.35);
	   int b = (int)(g * 0.35);
	   c.triangle(x, y, new Pixel(r, g, b));
	   }
	   }
	   }
	   c.save("out.ppm");
	   // ==================================== */
    }
}
