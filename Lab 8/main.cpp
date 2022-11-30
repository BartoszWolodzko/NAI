// opencv
#include "opencv2/opencv.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
// C++
#include <iostream>

using namespace cv;
using namespace std;

int main(int argc, char** argv) {
    VideoCapture cap(0);
    if (cap.isOpened() == false)
    {
        cout << "Cannot open the web cam" << endl;
        cin.get();
        return -1;
    }

    double dWidth = cap.get(CAP_PROP_FRAME_WIDTH);
    double dHeight = cap.get(CAP_PROP_FRAME_HEIGHT);

    cout << "Resolution of the video : " << dWidth << " x " << dHeight << endl;

    String windowNameOriginalImage = "Original Image";

    namedWindow(windowNameOriginalImage, WINDOW_NORMAL);
    namedWindow("Mirrored frame", WINDOW_NORMAL);

    while (true)
    {
        Mat frame;
        bool bSuccess = cap.read(frame);
        Mat frameMirror;
        flip(frame, frameMirror, 1);

        if (bSuccess == false) {
            cout << "Found the end of the video" << endl;
            break;
        }

        imshow(windowNameOriginalImage, frame);
        imshow("Mirrored frame", frameMirror);

        if (waitKey(10) == 27) {
            cout << "Esc key is pressed by user. Stoppig the video" << endl;
            break;
        }
    }
    return 0;
}