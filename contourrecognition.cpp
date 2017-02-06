#include "contourrecognition.h"

ContourRecognition::ContourRecognition(QObject *parent) : QObject(parent)
{

}

ContourRecognition::~ContourRecognition()
{

}

QImage ContourRecognition::getContour(QString imagePath)
{
    if(!imagePath.isEmpty() && !imagePath.isNull()){
        QImage image(imagePath);
        QColor color;
        QImage outImage(image.size(), QImage::Format_ARGB32);
        outImage.fill(Qt::white);

        int wMin = 50;
        int wMax = 280;
        int hMin = 230;
        int hMax = 340;

        for(int i = 0; i < image.width(); i++){
            for(int j = 0; j < image.height(); j++){
                color = image.pixelColor(i, j);
                if(color.blackF() > 0.3 && i > wMin && i < wMax && j > hMin && j < hMax)
                    outImage.setPixelColor(i, j, color);
            }
        }

        return outImage;
    }
    else
        return QImage(NULL);
}

QImage ContourRecognition::blowUpContour(QImage &image)
{
    if(!image.isNull()){
        QImage outImage = image;
        QPointF weightPoint(0,0);
        double num;

        for(int i = 0; i < image.width(); i++){
            for(int j = 0; j < image.height(); j++){
                if(image.pixelColor(i,j) != Qt::white){
                    weightPoint.setX(weightPoint.x() + i);
                    weightPoint.setY(weightPoint.y() + j);
                    num++;
                }
            }
        }
        weightPoint.setX(weightPoint.x() / num);
        weightPoint.setY(weightPoint.y() / num);

        for(double k = 0; k < image.width(); k++){
            for(double l = 0; l < image.height(); l++){
                if(image.pixelColor((int)k,(int)l) != Qt::white){
                    QPointF currentBlown = blownPoint(k-weightPoint.x(), l-weightPoint.y(), 10);
                    outImage.setPixelColor(qRound(currentBlown.x() + weightPoint.x()), qRound(currentBlown.y() + weightPoint.y()), Qt::magenta);
                }
            }
        }

        outImage.setPixelColor(qRound(weightPoint.x()), qRound(weightPoint.y()), Qt::yellow);
        return outImage;
    }
    else
        return QImage(NULL);
}

QPointF ContourRecognition::blownPoint(double x, double y, double blow)
{
    QPointF blown;
    blown.setX((sqrt((x*x+y*y)) + blow)*(x/(sqrt((x*x+y*y)))));
    blown.setY((sqrt((x*x+y*y)) + blow)*(y/(sqrt((x*x+y*y)))));
    return blown;
}

