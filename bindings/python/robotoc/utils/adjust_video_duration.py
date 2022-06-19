import sys
import os
import importlib


def get_video_duration_sec(video):
    cv2_spec = importlib.util.find_spec("cv2")
    found = cv2_spec is not None
    if found:
        import cv2
        cap = cv2.VideoCapture(video)
        return float(cap.get(cv2.CAP_PROP_FRAME_COUNT) / cap.get(cv2.CAP_PROP_FPS))
    else:
        return None


def adjust_video_duration(video: str, desired_duration_sec: float, output_name=None):
    video_rutation_sec = get_video_duration_sec(video)
    if video_rutation_sec is not None:
        video_name_split = video.rsplit('.', maxsplit=1)
        if output_name is None:
            output_name = video_name_split[0] + '_adjusted.' + video_name_split[1]
        ptsrate = video_rutation_sec / desired_duration_sec
        os.system('ffmpeg -i ' + video + ' -vf setpts=PTS/' + str(ptsrate) + ' -an ' + output_name)
    else:
        print('Could not find opencv-python!')


if __name__ == '__main__':
    video = sys.argv[1]
    desired_duration_sec = sys.argv[2]
    adjust_video_duration(video, desired_duration_sec)