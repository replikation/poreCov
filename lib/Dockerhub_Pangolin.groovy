import java.net.URL
import java.net.URLConnection

class Dockerhub_Pangolin {
    static boolean IsPangoAvailable() {
        try {
            final URL url = new URL("https://registry.hub.docker.com/v2/repositories/nanozoo/pangolin-v4/tags/");
            final URLConnection conn = url.openConnection();
            conn.connect();
            conn.getInputStream().close();
            return true;
        } catch (MalformedURLException e) {
            return false;
        } catch (IOException e) {
            return false;
        }
    }
}