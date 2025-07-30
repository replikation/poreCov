import java.net.URL
import java.net.URLConnection

class Git_PoreCov {
    
    static boolean IsGitAvailable() {
        try {
            final URL url = new URL("https://api.github.com/repos/replikation/poreCov/releases/latest");
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