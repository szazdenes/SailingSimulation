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
        int hMin = 280;
        int hMax = 340;

        for(int i = 0; i < image.width(); i++){
            for(int j = 0; j < image.height(); j++){
                color = image.pixelColor(i, j);
                if(color.blackF() > 0.3 && i > wMin && i < wMax && j > hMin && j < hMax)
                    outImage.setPixelColor(i, j, color);
            }
        }

        for(int k = 1; k < outImage.width()-1; k+=3){
            for(int l = 1; l < outImage.height()-1; l+=3){
                skeletonize(k, l, outImage);
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

        QList<QPoint> neighbour = getNeighbourList(image);

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

void ContourRecognition::skeletonize(int posX, int posY, QImage &image)
{
    double num = 0;
    double xPos = 0, yPos = 0;

    for(int i = posX-1; i <= posX+1; i++){
        for(int j = posY-1; j <= posY+1; j++){
            if(image.pixelColor(i,j) != Qt::white){
                xPos += i;
                yPos += j;
                num++;
            }
            image.setPixelColor(i,j, Qt::white);
        }
    }

    if(num != 0){
        xPos /= num;
        yPos /= num;
        QPoint result(qRound(xPos), qRound(yPos));
        image.setPixelColor(result, Qt::black);
    }
}

QList<QPoint> ContourRecognition::getNeighbourList(QImage &image)
{
    QList<QVector2D> pixelVector;

    for(int i = 0; i < image.height(); i++){
        for(int j = 0; j < image.width(); j++){
            if(image.pixelColor(j,i) != Qt::white)
                pixelVector.append(QVector2D(j,i));
        }
    }

    QMap<int, double> minDistanceMap;
    QList<QVector2D> neighbourVector;

    neighbourVector.append(pixelVector.at(0));

    for(int l = 0; l < pixelVector.size(); l++){
        for(int k = l; k < pixelVector.size(); k++){
            double dist = (pixelVector.at(k) - neighbourVector.last()).length();
            if(dist != 0 /*&& dist < 10*/)
                minDistanceMap[k] = dist;
        }

        QList<double> distanceList;

        foreach(double current, minDistanceMap.values())
            distanceList.append(current);

        qSort(distanceList);

        foreach(double current, distanceList){
            QVector2D currentVector = pixelVector.at(minDistanceMap.key(current));
            if(!neighbourVector.contains(currentVector)){
                neighbourVector.append(currentVector);
                break;
            }
        }
    }

    QList<QPoint> neighbourList;

    QImage neighbourImage(image.size(), QImage::Format_ARGB32);
    neighbourImage.fill(Qt::white);

    QPainter painter(&neighbourImage);
    painter.setPen(Qt::red);

    foreach(QVector2D current, neighbourVector)
        neighbourList.append(QPoint(qRound(current.x()), qRound(current.y())));


    for(int i = 0; i < neighbourList.size()-1; i++){
        painter.drawLine(neighbourList.at(i), neighbourList.at(i+1));
        neighbourImage.setPixelColor(neighbourList.at(i), Qt::black);

    }

    neighbourImage.save("../neighbour.png");

    return neighbourList;


}

QList<QPointF> ContourRecognition::scaleDataToImage(QString dataPath, QImage &image)
{
    QList<QPointF> dataList;
    QList<double> xList, yList;
    QFile file(dataPath);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        qDebug("Opening error.");
    }
    double x, y, xmin, xmax, ymin, ymax;
    QTextStream stream(&file);
    while(!stream.atEnd()){
        QString line = stream.readLine();
        QTextStream linestream(&line);
        linestream >> x >> y;
        dataList.append(QPointF(x, y));
        xList.append(x);
        yList.append(y);
    }

    file.close();

    qSort(xList);
    qSort(yList);
    xmin = xList.first();
    xmax = xList.last();
    ymin = yList.first();
    ymax = yList.last();

    QList<QPointF> scaledList;

    foreach(QPointF current, dataList){
        QPointF scaled;
        scaled.setX((image.width()/qAbs(xmax-xmin))*(current.x() - xmin));
        scaled.setY(-1*(image.height()/qAbs(ymax-ymin))*(current.y() - ymax));
        scaledList.append(scaled);
    }

//    qDebug("%f %f %f %f", xmax-xmin, ymax-ymin, xmin, ymax);
    return scaledList;
}

double ContourRecognition::blowDistance(double R, double s, double H, double h)
{
    double distance;
    double phi, theta;

    phi = acos(R/(R+H));
    theta = acos(R/(R+h));

    distance = R*(phi+theta) - s;

    return distance;
}

