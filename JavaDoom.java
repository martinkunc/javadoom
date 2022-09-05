import org.lwjgl.*;
import org.lwjgl.glfw.*;
import org.lwjgl.opengl.*;
import org.lwjgl.stb.STBEasyFont;
import org.lwjgl.system.MemoryUtil;


import java.nio.*;
import java.util.ArrayList;
import java.util.Collections;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import static org.lwjgl.glfw.Callbacks.*;
import static org.lwjgl.glfw.GLFW.*;
import static org.lwjgl.opengl.GL11.*;
import static org.lwjgl.system.MemoryUtil.*;

import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.Optional;

public class JavaDoom {

	// The window handle
	private static long window;
	private static FrameRate frameRate;

	/* Define window size */
	static int W = 608;
	static int H = 480;
	static float EyeHeight = 6f;  // Camera height from floor when standing
	static float DuckHeight = 2.5f;  // And when crouching
	static int HeadMargin = 1;    // How much room there is above camera before the head hits the ceiling
	static int KneeHeight = 2;    // How tall obstacles the player can simply walk over without jumping
	static double hfov = (0.73 * H);  // Affects the horizontal field of vision
	static double vfov = (.2 * H);    // Affects the vertical field of vision
	static boolean exit;
	static int[] wsad;
	static int ground;
    static int falling;
    static int moving;
    static int ducking;
    static float yaw;
	static double oldx = -1;
	static double oldy = -1;
	static boolean insideMouseCallback = false;
	static long lastDurationTime;

	 static class xy
    {
        public xy(Float float1, Float float2) {
			x = float1;
			y = float2;
		}

		public float x, y;
    }
     static class xyz
    {
		public xyz(float x, float y, float z) {
			this.x = x;
			this.y = y;
			this.z = z;
		}
        public float x, y, z;
    }
    /* Sectors: Floor and ceiling height; list of edge vertices and neighbors */
    static class sector
    {
        public float floor, ceil;
        public xy[] vertex; // Each vertex has an x and y coordinate
        public int[] neighbors; // Each edge may have a corresponding neighboring sector
        public int npoints;  // How many vertexes there are
    }
	static ArrayList<sector> sectors;

	static int NumSectors = 0;
    /* Player: location */
    static class Tplayer
    {
		public Tplayer(	xyz where, xyz velocity, float angle, int sector) 
		{
			this.where = where;
			this.velocity = velocity;
			this.angle = angle;
			this.sector = sector;
		}

        public xyz where;  // Current position,
		public xyz velocity;   // Current motion vector
        public float angle, anglesin, anglecos, yaw; // Looking towards (and sin() and cos() thereof)
        public int sector; // Which sector the player is currently in
    }
    static Tplayer player;

	// Clamps the value a among mi and ma
    static int clamp(int a, int mi, int ma)
    {
        return Math.min(Math.max(a, mi), ma);
    }

	static double clamp(double a, double mi, double ma)
	{
		return Math.min(Math.max(a, mi), ma);
	}

	// Returns cross product of vectors (x0,y0) and (x1, y1)
    static float vxs(float x0, float y0, float x1, float y1)
    {
        // vxs: Vector cross product
        return x0 * y1 - x1 * y0;
    }

	// Determines if two ranges a0-a1 and b0-b1 overlaps
    static boolean Overlap(float a0, float a1, float b0, float b1)
    {
        // Overlap:  Determine whether the two number ranges overlap.
        return (Math.min(a0, a1) <= Math.max(b0, b1) && Math.min(b0, b1) <= Math.max(a0, a1));
    }

	// Determines of two 2d boxes (x0, y0, x1, y1) and (x2, y2, x3, y3) intersect
    static boolean IntersectBox(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3)
    {
        return Overlap(x0, x1, x2, x3) && Overlap(y0, y1, y2, y3);
    }

	// Determines on which side of line of a line the points is on.
	// Returns value <0, =0 and >0
	// In pure math, cross product of 2d vector doesn't make sense, but is is extension of 
	// 3d vectors cross product where z coordinate is 0 (like 2d vectors in 3d spece with z axis=0).
	// It results in a vector which has only z coordinate.
	// It results in winding of vectors. If we go from a to b to c, will we be going clock wise or anti-clockwise.
	static float PointSide(float px, float py, float x0, float y0, float x1, float y1)
    {
        return vxs(x1 - x0, y1 - y0, px - x0, py - y0);
    }
	
	// Calculates the point of intersection between two lines.
    static xy Intersect(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4)
    {
		if (x1 == x2 && y1 == y2) return new xy(x1, y1);
        return new JavaDoom.xy(
			vxs(vxs(x1, y1, x2, y2), (x1) - (x2), vxs(x3, y3, x4, y4), (x3) - (x4)) / vxs((x1) - (x2), (y1) - (y2), (x3) - (x4), (y3) - (y4)),
			vxs(vxs(x1, y1, x2, y2), (y1) - (y2), vxs(x3, y3, x4, y4), (y3) - (y4)) / vxs((x1) - (x2), (y1) - (y2), (x3) - (x4), (y3) - (y4))
		);
        
    }

	private static Optional<Integer> IntTryParse(String string) {
		Optional<Integer> f = Optional.empty();
		try { f = Optional.of(Integer.parseInt(string)); } catch (Exception e) {}
		return f;
	}

	static void LoadData()
	{
		Pattern r = Pattern.compile("\\s+");
		try (BufferedReader sr = new BufferedReader(new FileReader("map-clear.txt"))) {
			String Buf;
			var NumVertices = 0;
			var vert = new ArrayList<xy>();
			sectors = new ArrayList<sector>();
			while ((Buf = sr.readLine()) != null)
			{
				var ba = r.split(Buf);
				if (ba.length == 0)
					continue;

				if (ba[0].equals("vertex"))
				{
					var vy = FloatTryParse(ba[1]);
					for (int i = 2; i <= ba.length - 1; i++)
					{
						var vx = FloatTryParse(ba[i]);
						if ( vy.isPresent() && vx.isPresent())
						{
							vert.add(new xy(vx.get(),vy.get() ));
							System.out.printf("Vertice %d (%f,%f) %n", NumVertices, vx.get(), vy.get());
							NumVertices++;
						}
					}
				}

				if (ba[0].equals("sector"))
				{
					NumSectors++;
					var sector = new sector();

					int i = 1;
					sector.floor = FloatTryParse(ba[i++]).get();
					sector.ceil = FloatTryParse(ba[i++]).get();
					var m = ba.length - i;

					sector.npoints = (int)m / 2;
					sector.vertex = new xy[m / 2 + 1];
					sector.neighbors = new int[m / 2];
					for (var j = 0; j < m / 2; j++)
					{
						var t = IntTryParse(ba[i + j]);
						if (t.isPresent())
						{
							sector.vertex[j + 1] = vert.get(t.get());
						}
					}
					sector.vertex[0] = sector.vertex[m/2];

					for (var j = 0; j < m / 2; j++)
					{
						var t = IntTryParse(ba[i + m / 2 + j]);
						if (t.isPresent())
						{
							sector.neighbors[j] = t.get();
						}
					}
					String neighborsS = IntStream.of(sector.neighbors).mapToObj(n -> String.valueOf(n)).collect(Collectors.joining(",")).toString();
					String verticesS = Stream.of(sector.vertex).map(n -> String.format("(%f, %f)", n.x, n.y)).collect(Collectors.joining(",")).toString();

					System.out.printf("Sector %d " +
									"Neighbours %s " +
									"Vertices %s %n",
							NumSectors - 1,
							neighborsS,
							verticesS);
					sectors.add(sector);
				}

				if (ba[0].equals("player"))
				{
					var vx = FloatTryParse(ba[1]);
					var vy = FloatTryParse(ba[2]);
					var angle = FloatTryParse(ba[3]);
					var sector = IntTryParse(ba[4]);
					if (vx.isPresent() &&
							vy.isPresent() &&
							angle.isPresent() &&
							sector.isPresent())
					{
						player = new Tplayer(
								new xyz( vx.get(), vy.get(),  0 ),
								new xyz( 0,  0, 0 ),
								angle.get(),
								(int)sector.get()
						);
						player.where.z = sectors.get((int)player.sector).floor + EyeHeight;
						System.out.printf("Player (%f,%f,%f) Angle %f, Sector %d %n", player.where.x, player.where.y, player.where.z, player.angle, player.sector);
					}
				}
			}
		} catch(IOException e) {}
	}

	private static Optional<Float> FloatTryParse(String string) {
		Optional<Float> f = Optional.empty();
		try { f = Optional.of(Float.parseFloat(string)); } catch (Exception e) {}
		return f;
	}
	static class Color {
		public Color(float r, float g,float b,float a) {
			R = r/(float)255;
			G = g/(float)255;
			B = b/(float)255;
			A = a/(float)255;
		}
		public float R;
        public float G;
        public float B;
        public float A;
		FloatBuffer buffer;
		public FloatBuffer toFloatButter() {
			var d = new float[]{R, G, B, A};
			buffer = MemoryUtil.memAllocFloat(d.length);
			buffer.put(d).flip();
			return buffer;
			//var b = BufferUtils.createFloatBuffer(4);
//			b.put(new float[]{R, G, B, A});
//			b.rewind();
//			return b;
		}

		public void free() {
			MemoryUtil.memFree(buffer);
		}
	}

	// Renders a vertical line at column x from y1 to y2 with color at top, middle and bottom
	static void vline(int ix, int iy1, int iy2, Color top, Color middle, Color bottom)
	{
		glLineWidth(1);
		var x = (double)ix;
		var y1 = (double)iy1;
		var y2 = (double)iy2;

		y1 = clamp(y1, 0, H - 1);
		y2 = clamp(y2, 0, H - 1);
		if (y2 == y1)
		{
			glBegin(GL_POINTS);
			glColor4fv(middle.toFloatButter());
			glVertex2d((double)x, (double)y1);
			glEnd();
		} 
		else if (y2 > y1)
		{
			glBegin(GL_POINTS);
			glColor4fv(top.toFloatButter());
			glVertex2d((double)x, (double)y1);
			glEnd();

//			glBegin(GL_POINTS);
//			glColor4fv(middle.toFloatButter());
//			for(int y=y1+1; y<y2; ++y) {
//				glVertex2d((double)x, (double)y);
//			}
//			glEnd();

			glBegin(GL_LINES);
			glColor4fv(middle.toFloatButter());
			glVertex2d((double)x, (double)y1 + 1);
			glVertex2d((double)x, (double)y2 - 1);
			glEnd();

			glBegin(GL_POINTS);
			glColor4fv(bottom.toFloatButter());
			glVertex2d((double)x, (double)y2);
			glEnd();

		}
	}

	// Renders a vertical line at column x from y1 to y2 with color at top, middle and bottom
	static void drawLine(int x1, int y1, int x2, int y2, Color color, int size)
	{
		glLineWidth(size);
		glBegin(GL_LINES);
		glColor4fv(color.toFloatButter());
		glVertex2d((double)x1, (double)y1);
		glVertex2d((double)x2, (double)y2);
		glEnd();
	}

	static void drawText(int x, int y, String t, Color color) {
		t = t + " ";
		ByteBuffer charBuffer = BufferUtils.createByteBuffer(t.length() * 270);
		int quads = STBEasyFont.stb_easy_font_print(x, y, t, null, charBuffer);

		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(2, GL_FLOAT, 16, charBuffer);
		glDrawArrays(GL_QUADS, 0, quads * 4);

	}


	/* MovePlayer(dx,dy): Moves the player by (dx,dy) in the map, and
     * also updates their anglesin/anglecos/sector properties properly.
     */
    static void MovePlayer(float dx, float dy)
    {
        float px = player.where.x; float py = player.where.y;
		// not used
//		var dx = 0;
//		var dy = 0;
        /* Check if this movement crosses one of this sector's edges
         * that have a neighboring sector on the other side.
         * Because the edge vertices of each sector are defined in
         * clockwise order, PointSide will always return -1 for a point
         * that is outside the sector and 0 or 1 for a point that is inside.
         */
        var sector = sectors.get((int)player.sector);
        var vert = sector.vertex;
        for(var s = 0; s< sector.npoints; ++s)
        {
            if (sector.neighbors[s] >=0 
                && IntersectBox(px, py, px + dx, py + dy, vert[s + 0].x, vert[s + 0].y, vert[s + 1].x, vert[s + 1].y)
                && PointSide(px + dx, py + dy, vert[s + 0].x, vert[s + 0].y, vert[s + 1].x, vert[s + 1].y) < 0)
            {
                player.sector = (int)sector.neighbors[s];
                break;
            }
        }

        player.where.x += dx;
        player.where.y += dy;
        player.anglesin = (float)Math.sin(player.angle);
        player.anglecos = (float)Math.cos(player.angle);
    }

	static class item {
		public item(int sectorNo, int sx1, int sx2) {
			this.sectorno = sectorNo;
			this.sx1 = sx1;
			this.sx2 = sx2;
		}
        public int sectorno;
        public int sx1;
        public int sx2;
    }

	static float Yaw(float y, float z, float yaw)
    {
        return y + z * yaw;
    }

	/*
 * initialize a smaller piece of the array and use the System.arraycopy 
 * call to fill in the rest of the array in an expanding binary fashion
 */
public static void arrayfill(int[] array, int value, int to) {
	int len = Math.min(array.length, to);
  
	//Value of i will be [1, 2, 4, 8, 16, 32, ..., len]
	for (int i = 0; i < len; i++) {
		array[i] = value;
	}
  }

	static void DrawScreen()
    {
// 32
        int MaxQueue = 16;  // maximum number of pending portal renders
        var queue = new item[MaxQueue];
        int head_pos = 0;
        int tail_pos = 0;
        var head = queue[head_pos];
        var tail = queue[tail_pos];
        var ytop = new int[W];
		arrayfill(ytop, 0, W);
        
        var ybottom = new int[W];
		// do we need this filling ?
        var renderedsectors = new ArrayList<Integer>( Collections.nCopies((int)NumSectors, 0));

        for (var x = 0; x < W; ++x) ybottom[x] = H - 1;
        for (var n = 0; n < NumSectors; ++n) renderedsectors.set(n, 0);

        /* Begin whole-screen rendering from where the player is. */
        queue[head_pos] = new item((int)player.sector,  0,  W - 1 );

        if (++head_pos == MaxQueue) head_pos = 0;

		long startTime = System.nanoTime();
        do
        {
            /* Pick a sector & slice from the queue to draw */
            item now = queue[tail_pos];
            if (++tail_pos == MaxQueue) tail_pos = 0;
			//System.out.printf("head_pos %d tail_pos %d  %n", tail_pos, head_pos);
			// only render 32 sectors per sector in queue, then jumps to next item from queue
            if ((renderedsectors.get(now.sectorno) & 0x21) > 0) continue; // Odd = still rendering, 0x20 = give up
            renderedsectors.set(now.sectorno, renderedsectors.get(now.sectorno)+1);
            sector sect = sectors.get(now.sectorno);

            /* Render each wall of this sector that is facing towards player. */
            // npoints is number of neighbour vectors in current sector
            for(int s=0;s < sect.npoints; ++s)
            {
                /* Acquire the x,y coordinates of the two endpoints (vertices) of this edge of the sector */
				/* Move the vector by player's position */
                float vx1 = sect.vertex[s + 0].x - player.where.x, vy1 = sect.vertex[s + 0].y - player.where.y;
                float vx2 = sect.vertex[s + 1].x - player.where.x, vy2 = sect.vertex[s + 1].y - player.where.y;
                /* Rotate them around the player's view */
				/* Rotates vectors vx and vy by angle player.angle looking from top  */
				float pcos = player.anglecos, psin = player.anglesin;
				// vx, vy -> tx, tz
                float tx1 = vx1 * psin - vy1 * pcos, tz1 = vx1 * pcos + vy1 * psin;
                float tx2 = vx2 * psin - vy2 * pcos, tz2 = vx2 * pcos + vy2 * psin;
                /* Is the wall at least partially in front of the player? */
				// = =
                if (tz1 < 0 && tz2 < 0) continue;

                /* If it's partially behind the player, clip it against player's view frustum */
                if (tz1 <= 0 || tz2 <= 0)
                {
                    float nearz = 1e-4f, farz = 5, nearside = 1e-5f, farside = 20f;
                    // Find an intersection between the wall and the approximate edges of player's view
                    xy i1 = Intersect(tx1, tz1, tx2, tz2, -nearside, nearz, -farside, farz);
                    xy i2 = Intersect(tx1, tz1, tx2, tz2, nearside, nearz, farside, farz);
                    
                    if (tz1<nearz) { if (i1.y > 0) { tx1 = i1.x; tz1 = i1.y; } else { tx1 = i2.x; tz1 = i2.y; } }
                    if (tz2<nearz) { if (i1.y > 0) { tx2 = i1.x; tz2 = i1.y; } else { tx2 = i2.x; tz2 = i2.y; } }
                }
                /* Do perspective transformation */
                float xscale1 = (float)hfov / tz1, yscale1 = (float)vfov / tz1; int x1 = W / 2 - (int)(tx1 * xscale1);
                float xscale2 = (float)hfov / tz2, yscale2 = (float)vfov / tz2; int x2 = W / 2 - (int)(tx2 * xscale2);
                if (x1 >= x2 || x2 < now.sx1 || x1 > now.sx2) continue; // Only render if it's visible
                /* Acquire the floor and ceiling heights, relative to where the player's view is */
                float yceil = sect.ceil - player.where.z;
                float yfloor = sect.floor - player.where.z;
                /* Check the edge type. neighbor=-1 means wall, other=boundary between two sectors. */
                int neighbor = sect.neighbors[s];
                float nyceil = 0, nyfloor = 0;
                if (neighbor >= 0) // Is another sector showing through this portal?
                {
                    nyceil = sectors.get(neighbor).ceil - player.where.z;
                    nyfloor = sectors.get(neighbor).floor - player.where.z;
                }
                /* Project our ceiling & floor heights into screen coordinates (Y coordinate) */
                int y1a = H / 2 - (int)(Yaw(yceil, tz1, player.yaw) * yscale1), y1b = H / 2 - (int)(Yaw(yfloor, tz1, player.yaw) * yscale1);
                int y2a = H / 2 - (int)(Yaw(yceil, tz2, player.yaw) * yscale2), y2b = H / 2 - (int)(Yaw(yfloor, tz2, player.yaw) * yscale2);
                /* The same for the neighboring sector */
                int ny1a = H / 2 - (int)(Yaw(nyceil, tz1, player.yaw) * yscale1), ny1b = H / 2 - (int)(Yaw(nyfloor, tz1, player.yaw) * yscale1);
                int ny2a = H / 2 - (int)(Yaw(nyceil, tz2, player.yaw) * yscale2), ny2b = H / 2 - (int)(Yaw(nyfloor, tz2, player.yaw) * yscale2);

                /* Render the wall. */
                int beginx = Math.max(x1, now.sx1), endx = Math.min(x2, now.sx2);
                for(int x = beginx; x <= endx; ++x)
                {
                    /* Calculate the Z coordinate for this point. (Only used for lighting.) */
                    int z = (int)(((x - x1) * (tz2 - tz1) / (x2 - x1) + tz1) * 8);

                    /* Acquire the Y coordinates for our ceiling & floor for this X coordinate. Clamp them. */
                    int ya = (x - x1) * (y2a - y1a) / (x2 - x1) + y1a, cya = clamp(ya, ytop[x], ybottom[x]); // top
                    int yb = (x - x1) * (y2b - y1b) / (x2 - x1) + y1b, cyb = clamp(yb, ytop[x], ybottom[x]); // bottom

                    /* Render ceiling: everything above this sector's ceiling height. */
					var c1Color = new Color(17,17,17,255); /* 0x111111 */
					var c2Color = new Color(54,54,54,255); /* 0x222222 */
					var c3Color = new Color(17,17,17,255); /* 0x111111 */
					vline(x, ytop[x], cya - 1, c1Color, c2Color, c3Color) ;
					c1Color.free();
					c2Color.free();
					c3Color.free();
            
                    /* Render floor: everything below this sector's floor height. */
					var f1Color = /* 0x0000FF */ new Color(0,0,255,255);
					var f2Color = /* 0x0000AA */ new Color(0,0,170,255);
					var f3Color = /* 0x0000FF */ new Color(0,0,255,255);
					vline(x, cyb+1, ybottom[x],f1Color,f2Color,f3Color);
					f1Color.free();
					f2Color.free();
					f3Color.free();

                    /* Is there another sector behind this edge? */
                    if (neighbor >= 0)
                    {
                        /* Same for _their_ floor and ceiling */
                        int nya = (x - x1) * (ny2a - ny1a) / (x2 - x1) + ny1a, cnya = clamp(nya, ytop[x], ybottom[x]);
                        int nyb = (x - x1) * (ny2b - ny1b) / (x2 - x1) + ny1b, cnyb = clamp(nyb, ytop[x], ybottom[x]);
                        /* If our ceiling is higher than their ceiling, render upper wall */
                        var r1 = new Color(1 * 255 -z,1 * 255 -z,1 * 255 -z, 255);
						//var r1 = new Color(0,1 * 255 -z,0, 255);
                        // 0x040007 * (31-z/8);
                        var r2 = new Color(4 * 31 - z / 8,0,7 * 31 - z / 8, 255);
						// Ceiling to neighbour ceiling line
                        vline(x, cya,cnya - 1, new Color(0,0,0, 255), x == x1 || x == x2 ? new Color( 0,0,0, 255) : r1, new Color( 0,0,0, 255)); // Between our and their ceiling
                        ytop[x] = clamp(Math.max(cya, cnya), ytop[x], H - 1);   // Shrink the remaining window below these ceilings
                        /* If our floor is lower than their floor, render bottom wall */
                        vline(x, cnyb + 1, cyb, new Color(0,0,0,255), x == x1 || x == x2 ? new Color( 0,0,0, 255) : r2, new Color(0,0,0, 255)); // Between their and our floor
                        ybottom[x] = clamp(Math.min(cyb, cnyb), 0, ybottom[x]); // Shrink the remaining window above these floors
                    }
                    else
                    {
                        /* There's no neighbor. Render wall from top (cya = ceiling level) to bottom (cyb = floor level). */
                        //unsigned r = 0x010101 * (255-z);
                        Color r = new Color( 1 * 255 - z, 1 * 255 - z, 1 * 255 - z, 255);
                        vline(x, cya, cyb, new Color( 0,0,0, 255), x == x1 || x == x2 ? new Color(0,0,0,255) : r, new Color( 0,0,0, 255));
                    }
                } // end for of wall rendering

                /* Schedule the neighboring sector for rendering within the window formed by this wall. */
                if (neighbor >= 0 && endx >= beginx && ((head_pos + MaxQueue + 1 - tail_pos) % MaxQueue > 0))
                {
                    queue[head_pos] = new item( neighbor, beginx, endx);
                                    
                    if (++head_pos == MaxQueue) head_pos = 0;
                }
            } // end of for for sector

            renderedsectors.set(now.sectorno, renderedsectors.get(now.sectorno)+1);
        } while (head_pos!=tail_pos);
		lastDurationTime = System.nanoTime() - startTime;
		//DrawMap(player);
    }

	public static void DrawMap(Tplayer player) {
		Color mapColor = new Color(255,255,255,255);
		Color playerColor = new Color(0,255,0,255);
		var mapTop = H / 2 - 200;
		var mapLeft = 10;
		int mapSquare = 20;
		var playerSize = 4;
		int wallThickness = 1;
		int playerThickness = 2;

		var secToRender = new ArrayList<Integer>();
		secToRender.add(player.sector);
		secToRender.addAll(IntStream.of(sectors.get(player.sector).neighbors).boxed().collect(Collectors.toList()));
		//Collections.addAll(secToRender, sectors.get(sector).neighbors);
		var playerSectorCeil = sectors.get(player.sector).ceil;
		for(var s = 0; s < sectors.size() ; s++) {
			//if (s < 0) continue;
			for(var i = 0; i < sectors.get(s).vertex.length; i++) {
				//
				//if (sectors.get(s).floor >  playerSectorCeil - EyeHeight) continue;

				var cv = sectors.get(s).vertex[i];
				xy nv;
				if (i < sectors.get(s).vertex.length - 1) {
					nv = sectors.get(s).vertex[i + 1];
					drawLine((int) cv.x * mapSquare + mapLeft, (int) cv.x * mapSquare + mapTop,
							(int) nv.x * mapSquare + mapLeft, (int) nv.y * mapSquare + mapTop, mapColor, wallThickness);
					drawText((int)nv.x * mapSquare + mapLeft, (int)nv.y* mapSquare + mapTop,
							String.format("%d %d %d", (int)nv.x, (int)nv.y, (int)sectors.get(s).floor), mapColor );
					drawText((int)cv.x * mapSquare + mapLeft, (int)cv.y* mapSquare + mapTop,
							String.format("%d %d %d", (int)cv.x, (int)cv.y, (int)sectors.get(s).floor), mapColor );
				}
//				else {
//					drawText((int)cv.x * mapSquare + mapLeft, (int)cv.y* mapSquare + mapTop, String.format("%d %d", (int)cv.x, (int)cv.y), mapColor );
//				}
			}
		}

		// player

		//drawArrow(player, playerColor, mapTop, mapLeft, mapSquare, playerSize, playerThickness);
		drawCross(player, playerColor, mapTop, mapLeft, mapSquare, playerSize, playerThickness);
	}

	private static void drawCross(Tplayer player, Color playerColor, int mapTop, int mapLeft, int mapSquare, int playerSize, int playerThickness) {
		drawLine((int) player.where.x* mapSquare - playerSize + mapLeft, (int) player.where.y* mapSquare - playerSize + mapTop,
				(int) player.where.x* mapSquare + playerSize + mapLeft, (int) player.where.y* mapSquare + playerSize + mapTop, playerColor, playerThickness);
		drawLine((int) player.where.x* mapSquare + playerSize + mapLeft, (int) player.where.y* mapSquare - playerSize + mapTop,
				(int) player.where.x* mapSquare - playerSize + mapLeft, (int) player.where.y* mapSquare + playerSize + mapTop, playerColor, playerThickness);
	}

	private static void drawArrow(Tplayer player, Color playerColor, int mapTop, int mapLeft, int mapSquare, int playerSize, int playerThickness) {
		var pAngle = player.angle;
		var angleSin = (float)Math.sin(-pAngle);
		var angleCos = (float)Math.cos(-pAngle);
		float px1 = (player.where.x - playerSize) * angleSin - player.where.y * angleCos;
		float py1 = (player.where.x - playerSize) * angleCos + player.where.y * angleSin;

		drawLine((int) player.where.x * mapSquare + mapLeft, (int) player.where.y * mapSquare + mapTop,
				(int) px1 * mapSquare + mapLeft, (int) py1 * mapSquare + mapTop, playerColor, playerThickness);
//		drawLine((int) player.where.x* mapSquare + playerSize + mapLeft, (int) player.where.y* mapSquare - playerSize + mapTop,
//				(int) player.where.x* mapSquare - playerSize + mapLeft, (int) player.where.y* mapSquare + playerSize + mapTop, playerColor, playerThickness);
	}

	private static void mouseCallback(long window, double x, double y) {
		if (insideMouseCallback) return;
		insideMouseCallback = true;
		IntBuffer leftBuff = BufferUtils.createIntBuffer(1);
		IntBuffer topBuff = BufferUtils.createIntBuffer(1);
		GLFW.glfwGetWindowPos(window, leftBuff, topBuff);

		int left = leftBuff.get(0);
		int top = topBuff.get(0);
		final IntBuffer bWidth = BufferUtils.createIntBuffer(1);
		final IntBuffer bHeight = BufferUtils.createIntBuffer(1);
		GLFW.glfwGetFramebufferSize(window, bWidth, bHeight);
		int width = bWidth.get(0);
		int height = bHeight.get(0);

//		var nx = (((left + x) / width) * 100) - 50;
//		var ny = (((top + y) / height) * 10) - 5;
		if (oldx == -1 && oldy == -1) {
			oldx = (int)x;
			oldy = (int)y;
		}
		var nx = clamp(((oldx - x) * 0.05) , -0.5,0.5);
		var ny = clamp(((oldy - y) * 0.05), -0.5 ,0.5);
		System.out.printf("nx %f, ny %f%n", nx, ny);
		oldx = x;
		oldy = y;
		// relative x mouse position * 0.05f
		player.angle = player.angle - (float)(nx);
		// yaw +
		//var de = clamp((int)(y * 0.05), -50, 50);
		//yaw = clamp((int)( yaw + y * 0.5), -10, 10);
		//player.yaw = yaw - player.velocity.z * 0.5f;
		player.yaw = player.yaw - (float)(ny);

		MovePlayer(0,0);

		//glfwSetCursorPos(window, 0, 0);
		insideMouseCallback = false;
	}
	public JavaDoom() {
		frameRate = new FrameRate();
	}
	public static void run() {

		try {
            
			// will print the error message in System.err.
			GLFWErrorCallback.createPrint(System.err).set();
	
			// Initialize GLFW. Most GLFW functions will not work before doing this.
			if ( !glfwInit() )
				throw new IllegalStateException("Unable to initialize GLFW");


			// Configure GLFW
			// glfwDefaultWindowHints(); // optional, the current window hints are already the default
			// glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE); // the window will stay hidden after creation
			// glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE); // the window will be resizable
			
			LoadData();
			
			wsad = new int[] { 0,0,0,0};
			 ground = 0;
             falling = 1;
             moving = 0;
             ducking = 0;
             yaw = 0;

			int w = W;
			int h = H;
			// Create the window
			window = glfwCreateWindow(w,h, "Java Doom", NULL, NULL);
			if ( window == NULL )
				throw new RuntimeException("Failed to create the GLFW window");

			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

			// Setup a key callback. It will be called every time a key is pressed, repeated or released.
			glfwSetKeyCallback(window, (window, key, scancode, action, mods) -> {
				if ( key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE ) {
					exit = true;
					glfwSetWindowShouldClose(window, true); // We will detect this in the rendering loop
				}
				if ( key == GLFW_KEY_W) {
					wsad[0] = action == GLFW_PRESS || action == GLFW_REPEAT ? 1 : 0;
				}
				if ( key == GLFW_KEY_S) {
					wsad[1] = action == GLFW_PRESS || action == GLFW_REPEAT ? 1 : 0;
				}
				if ( key == GLFW_KEY_A) {
					wsad[2] = action == GLFW_PRESS || action == GLFW_REPEAT ? 1 : 0;
				}
				if ( key == GLFW_KEY_D) {
					wsad[3] = action == GLFW_PRESS || action == GLFW_REPEAT ? 1 : 0;
				}
				if ( key == GLFW_KEY_SPACE) {
					if (ground == 1) { player.velocity.z += 0.5f; falling = 1; }
				}
				if ( key == GLFW_KEY_LEFT_CONTROL || key == GLFW_KEY_RIGHT_CONTROL) {
					ducking = action == GLFW_PRESS || action == GLFW_REPEAT ? 1 : 0;
                    falling = 1; 
				}

			});

			yaw = clamp((int)( yaw + 0 * 0.5), -10, 10);
			player.yaw = yaw - player.velocity.z * 0.5f;
			MovePlayer(0,0);

			glfwSetCursorPosCallback(window, (window, x, y) -> mouseCallback(window, x, y));

			
	
			// Get the thread stack and push a new frame
			// try ( MemoryStack stack = stackPush() ) {
			// 	IntBuffer pWidth = stack.mallocInt(1); // int*
			// 	IntBuffer pHeight = stack.mallocInt(1); // int*
	
			// 	// Get the window size passed to glfwCreateWindow
			// 	glfwGetWindowSize(window, pWidth, pHeight);
	
			// 	// Get the resolution of the primary monitor
			// 	GLFWVidMode vidmode = glfwGetVideoMode(glfwGetPrimaryMonitor());
	
			// 	// Center the window
			// 	glfwSetWindowPos(
			// 		window,
			// 		(vidmode.width() - pWidth.get(0)) / 2,
			// 		(vidmode.height() - pHeight.get(0)) / 2
			// 	);
			// } // the stack frame is popped automatically
	
			// Make the OpenGL context current
			glfwMakeContextCurrent(window);
			GL.createCapabilities();

			String glVersion = GL11.glGetString(GL11.GL_VERSION);
			System.out.println(glVersion);

			// glfwSetWindowRefreshCallback(window, windowHnd -> {
			// 	//gears.render();
			// 	//gears.animate();
			// 	glfwSwapBuffers(windowHnd);
			// });
			
			// Enable v-sync
			//glfwSwapInterval(1);
	
			// Make the window visible
			//glfwShowWindow(window);

	
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glOrtho( 0, w, h,0, -1, 1);
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			// https://www.gamedev.net/forums/topic/447084-glortho-off-by-one-errors/3958345/?topic_id=447084&whichpage=1
			glTranslatef(0.375f, 0.375f, 0);

			// // Initialization code OpenGL
			// glMatrixMode(GL_PROJECTION);
			// glLoadIdentity();
			// //glOrtho(0, 640, 480, 0, 1, -1);
			// glOrtho(-1, 1, -1, 1, -1, 1);
			// glMatrixMode(GL_MODELVIEW);
			
			// Disable vsync
			//glfwSwapInterval(0);
			
			
            

			frameRate.init();

			exit = false;
			while (!exit ) {
				if (glfwWindowShouldClose(window))
					exit = true;
				// Render
				frameRate.update(player, sectors, lastDurationTime);
				
				glClear(GL_COLOR_BUFFER_BIT);
				
				//GL11.glBindTexture(GL11.GL_TEXTURE_2D, tex);

				DrawScreen();
	
				glColor3f(1,1,1);
				frameRate.draw();
				

				float eyeheight = ducking == 1 ? DuckHeight : EyeHeight;
                ground = falling != 1 ? 1 : 0;
                if (falling==1)
                {
                    player.velocity.z -= 0.05f; // add gravity;
                    float nextz = player.where.z + player.velocity.z;
                    
                    // When going down
                    if (player.velocity.z < 0 && nextz < sectors.get((int)player.sector).floor + eyeheight)
                    {
                        /* Fix to ground */
                        player.where.z = sectors.get((int)player.sector).floor + eyeheight;
                        player.velocity.z = 0;
                        falling = 0;
                        ground = 1;
                    }
                    else if (player.velocity.z > 0 && nextz > sectors.get((int)player.sector).ceil) // When going up
                    {
                        /* Prevent jumping above ceiling */
                        player.velocity.z = 0;
                        falling = 1;
                    }
                    if (falling==1)
                    {
                        player.where.z += player.velocity.z;
                        moving = 1;
                    }
                }

                /* Horizontal collision detection */
                if (moving == 1)
                {
                    float px = player.where.x, py = player.where.y;
                    float dx = player.velocity.x, dy = player.velocity.y;

                    var sect = sectors.get((int)player.sector);
                    xy[] vert = sect.vertex;
                    /* Check if the player is about to cross one of the sector's edges */
                    for (var s = 0; s < sect.npoints; ++s)
                    {
                        if (IntersectBox(px, py, px + dx, py + dy, vert[s + 0].x, vert[s + 0].y, vert[s + 1].x, vert[s + 1].y)
                            && PointSide(px + dx, py + dy, vert[s + 0].x, vert[s + 0].y, vert[s + 1].x, vert[s + 1].y) < 0)
                        {
                            /* Check where the hole is. */
                            float hole_low = sect.neighbors[s] < 0 ? 9e9f : Math.max(sect.floor, sectors.get(sect.neighbors[s]).floor);
                            float hole_high = sect.neighbors[s] < 0 ? -9e9f : Math.min(sect.ceil, sectors.get(sect.neighbors[s]).ceil);
                            /* Check whether we're bumping into a wall. */
                            if (hole_high < player.where.z + HeadMargin
                            || hole_low > player.where.z - eyeheight + KneeHeight)
                            {
                                /* Bumps into a wall! Slide along the wall. */
                                /* This formula is from Wikipedia article "vector projection". */
                                float xd = vert[s + 1].x - vert[s + 0].x, yd = vert[s + 1].y - vert[s + 0].y;
                                dx = xd * (dx * xd + yd * dy) / (xd * xd + yd * yd);
                                dy = yd * (dx * xd + yd * dy) / (xd * xd + yd * yd);
                                moving = 0;
                            }
                        }
                    }

                    MovePlayer(dx,dy);
                    falling = 1;
                }

				float[] move_vec = new float[] { 0f, 0f };
                if (wsad[0] == 1) { move_vec[0] += player.anglecos * 0.2f; move_vec[1] += player.anglesin * 0.2f; }
                if (wsad[1] == 1) { move_vec[0] -= player.anglecos * 0.2f; move_vec[1] -= player.anglesin * 0.2f; }
                if (wsad[2] == 1) { move_vec[0] += player.anglesin * 0.2f; move_vec[1] -= player.anglecos * 0.2f; }
                if (wsad[3] == 1) { move_vec[0] -= player.anglesin * 0.2f; move_vec[1] += player.anglecos * 0.2f; }
                int pushing = wsad[0] == 1 || wsad[1] == 1 || wsad[2] == 1 || wsad[3] == 1 ? 1 : 0;
                float acceleration = pushing == 1 ? 0.4f : 0.2f;

                player.velocity.x = player.velocity.x * (1 - acceleration) + move_vec[0] * acceleration;
                player.velocity.y = player.velocity.y * (1 - acceleration) + move_vec[1] * acceleration;

                if (pushing == 1) moving = 1;
				
				glfwSwapBuffers(window);
				glfwPollEvents();
			}
	
			// Free the window callbacks and destroy the window
		glfwFreeCallbacks(window);
		glfwDestroyWindow(window);
		
		// Terminate GLFW and free the error callback
		glfwTerminate();
		glfwSetErrorCallback(null).free();
	} catch (Exception e) {
		e.printStackTrace();

		// Free the window callbacks and destroy the window
		glfwFreeCallbacks(window);
		glfwDestroyWindow(window);

		// Terminate GLFW and free the error callback
		glfwTerminate();
		glfwSetErrorCallback(null).free();
		
		System.exit(1);
	}
	}



}