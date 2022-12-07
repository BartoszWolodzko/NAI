#include <iostream>
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <vector>

using namespace cv;
using namespace std;

int main(){
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
    createTrackbar("Threshold1", windowNameOriginalImage, 0, 255, NULL);
    createTrackbar("Threshold2", windowNameOriginalImage, 0, 255, NULL);


    while (true)
    {
        Mat frame;
        bool bSuccess = cap.read(frame);

        if (!bSuccess) {
            cout << "Found the end of the video" << endl;
            break;
        }

        // cut all colors except red
        Mat redOnly;
        inRange(frame, Scalar(0, 0, 100), Scalar(100, 100, 255), redOnly);

        blur(redOnly, redOnly, Size(5, 5));

        Mat cannyImage;
        int threshold1 = getTrackbarPos("Threshold1", windowNameOriginalImage);
        int threshold2 = getTrackbarPos("Threshold2", windowNameOriginalImage);
        Canny(redOnly, cannyImage, threshold1, threshold2);

        // get contours from canny image
        vector<vector<Point>> contours;
        vector<Vec4i> hierarchy;
        findContours(cannyImage, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE, Point(0, 0));

        vector<Point> centers;
        for (const auto & contour : contours)
        {
            // calculate moments
            Moments moment = moments(contour);
            double area = moment.m00;
            if (area > 1000) {
                // calculate center of mass
                double x = moment.m10 / area;
                double y = moment.m01 / area;

                centers.emplace_back(x, y);
            }
        }

        for (int i = 0; i < centers.size(); i++)
        {
            for (int j = i + 1; j < centers.size(); j++)
            {
                // draw lines if centers are horizontally
                if (abs(centers[i].y - centers[j].y) < 10)
                {
                    // calculate line thickness based on distance between centers
                    int thickness = abs(centers[i].x - centers[j].x);
                    thickness=thickness/10;
                    cout << thickness<<endl;
                    if (thickness<=0) thickness=1;



                    line(frame, centers[i], centers[j], Scalar(0, 0, 255), 5*thickness);
                    line(frame, Point (centers[i].x, centers[i].y+5),
                         Point (centers[j].x, centers[j].y+5), Scalar(255, 0, 0), thickness);
                }
            }
        }


        imshow(windowNameOriginalImage, frame);

        if (waitKey(10) == 27) {
            cout << "Esc key is pressed by user. Stoppig the video" << endl;
            break;
        }
    }

    return 0;
}