import static org.lwjgl.glfw.GLFW.*;
import static org.lwjgl.opengl.GL11.*;

import java.nio.ByteBuffer;
import java.util.ArrayList;

import org.lwjgl.BufferUtils;
import org.lwjgl.stb.STBEasyFont;

public class FrameRate {
    double previousTime;
    int frameCount;
    String text;
    long lastDuration;
    final int DURATIONS_COUNT = 400;
    final int MAX_DURATION = 200;
    int[] durations = new int[DURATIONS_COUNT];
    int ctm;

    // flat whether to draw rendering routine duration in a line graph
    boolean drawHeartBeatEnabled;
    public void init() {
        previousTime = glfwGetTime();
		frameCount = 0;
		text = "Frames:";
        ctm = 0;
        drawHeartBeatEnabled = false;
    }

    public void update(JavaDoom.Tplayer player, ArrayList<JavaDoom.sector> sectors, long lastDurationTime) {
        double currentTime = glfwGetTime();
        lastDuration = lastDurationTime;
        frameCount++;
        // If a second has passed.
        if ( currentTime - previousTime >= 1.0 )
        {
            text = String.format("Frames: %d %n Player (%f,%f,%f) angle %f yaw %f%nsector %d floor %f ceil %f ", frameCount,
                    player.where.x, player.where.y, player.where.z, (player.angle / (2 * Math.PI)) * 360, player.yaw,
                    player.sector, sectors.get(player.sector).floor, sectors.get(player.sector).ceil);

            frameCount = 0;
            previousTime = currentTime;
        }
    }

    public void draw() {
        double currentTime = glfwGetTime();
        ByteBuffer charBuffer = BufferUtils.createByteBuffer(text.length() * 270);
        int quads = STBEasyFont.stb_easy_font_print(0, 0, text, null, charBuffer);

        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(2, GL_FLOAT, 16, charBuffer);
        glDrawArrays(GL_QUADS, 0, quads * 4);
        drawHeartBeat();

    }

    private void drawHeartBeat() {
        if (!drawHeartBeatEnabled)
            return;
        ctm = (int)(++ctm % DURATIONS_COUNT);
        durations[ctm] = (int)lastDuration;
        //if ( currentTime - previousTime >= 1.0 )
        for(int i=0; i < DURATIONS_COUNT; i++)
        {
            glBegin(GL_LINES);
            //glColor4fv(color.toFloatButter());
            glVertex2d((double)i , (double)0);
            glVertex2d((double)i, (double) Math.min(durations[i] / 1000000, MAX_DURATION));
            glEnd();
        }
    }
}
