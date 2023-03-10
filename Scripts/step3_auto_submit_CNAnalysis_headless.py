# coding=utf-8
from selenium import webdriver
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
from time import sleep
import re  # 用于正则
from PIL import Image  # 用于打开图片和对图片处理
import pytesseract  # 用于图片转文字
from selenium import webdriver  # 用于打开网站
import time  # 代码运行停顿
import os, sys


class VerificationCode:
    def __init__(self):
        # self.driver = webdriver.Chrome(executable_path="/usr/local/bin/chromedriver")
        self.driver = webdriver.PhantomJS(executable_path="/usr/local/bin/phantomjs")
        # self.driver = webdriver.Remote(desired_capabilities=DesiredCapabilities.HTMLUNIT)
        # self.driver = webdriver.DesiredCapabilities.HTMLUNIT
        self.find_element = self.driver.find_element_by_css_selector

    def get_pictures(self):
        self.driver.get('https://cpclab.uni-duesseldorf.de/cna/main.php')  # 打开登陆页面
        self.driver.save_screenshot('pictures.png')  # 全屏截图
        page_snap_obj = Image.open('pictures.png')
        img = self.find_element('#captcha')  # 验证码元素位置
        time.sleep(1)
        location = img.location
        size = img.size  # 获取验证码的大小参数
        left = location['x']
        top = location['y']
        right = left + size['width']
        bottom = top + size['height']
        image_obj = page_snap_obj.crop((left, top, right, bottom))  # 按照验证码的长宽，切割验证码
        # image_obj.show()  # 打开切割后的完整验证码
        # self.driver.close()  # 处理完验证码后关闭浏览器
        return image_obj

    def processing_image(self):
        image_obj = self.get_pictures()  # 获取验证码
        img = image_obj.convert("L")  # 转灰度
        pixdata = img.load()
        w, h = img.size
        threshold = 160
        # 遍历所有像素，大于阈值的为黑色
        for y in range(h):
            for x in range(w):
                if pixdata[x, y] < threshold:
                    pixdata[x, y] = 0
                else:
                    pixdata[x, y] = 255
        return img

    def delete_spot(self):
        images = self.processing_image()
        data = images.getdata()
        w, h = images.size
        black_point = 0
        for x in range(1, w - 1):
            for y in range(1, h - 1):
                mid_pixel = data[w * y + x]  # 中央像素点像素值
                if mid_pixel < 50:  # 找出上下左右四个方向像素点像素值
                    top_pixel = data[w * (y - 1) + x]
                    left_pixel = data[w * y + (x - 1)]
                    down_pixel = data[w * (y + 1) + x]
                    right_pixel = data[w * y + (x + 1)]
                    # 判断上下左右的黑色像素点总个数
                    if top_pixel < 10:
                        black_point += 1
                    if left_pixel < 10:
                        black_point += 1
                    if down_pixel < 10:
                        black_point += 1
                    if right_pixel < 10:
                        black_point += 1
                    if black_point < 1:
                        images.putpixel((x, y), 255)
                    black_point = 0
        #images.show()
        return images

    def image_str(self):
        image = self.delete_spot()
        pytesseract.pytesseract.tesseract_cmd = '/usr/bin/tesseract'  # 设置pyteseract路径r"C:\Program Files\Tesseract-OCR\tesseract.exe"

        result = pytesseract.image_to_string(image)  # 图片转文字
        print(result)
        resultj = re.sub(u"([^\u4e00-\u9fa5\u0030-\u0039\u0041-\u005a\u0061-\u007a])", "", result)  # 去除识别出来的特殊字符
        print(resultj)
        if len(resultj) > 4:
            result_four = resultj[1:5]  # 只获中间4个字符，会更适合这个网页
        if len(resultj) <= 4:
            result_four = resultj[0:4]  #取前四个字符
        # print(resultj)  # 打印识别的验证码
        if len(result_four) != 4:
            self.driver.quit()
        return result_four

    def dir_exists(self, dir):
        if os.path.exists(dir):
            print(dir)
        #        shutil.rmtree(dir)
        if not os.path.exists(dir):
            os.mkdir(dir)

    def get_web_data(self, submit_success_url):
        sleep(30)
        # self.driver = webdriver.Chrome(executable_path="/usr/local/bin/chromedriver")
        self.driver = webdriver.PhantomJS(executable_path="/usr/local/bin/phantomjs")
        self.find_element = self.driver.find_element_by_css_selector
        self.driver.get(submit_success_url)
        urls = self.driver.find_elements_by_xpath("//a")
        cna_results_dir = 'cna_results'
        self.dir_exists(cna_results_dir)
        for url in urls:
            url_link = url.get_attribute("outerHTML")
            # url_link = url.get_attribute("href").toString()
            if '.pdb' in url_link or '.dat' in url_link:
                download_link = url_link.split('>')[0].split('=')[-1].strip().strip('"')
                print(submit_success_url + '/' + download_link)
                os.system('wget -np -P ./' + cna_results_dir + ' ' + submit_success_url + '/' + download_link)
        self.driver.quit()

    def fill_data(self, result_four, analysis_type, pdbid, pdb_file, cut_hydro):
        # 查找页面的Aanlysis type里的“Single network”选项,他的id='m1'，并进行点击
        self.driver.find_elements_by_id(analysis_type)[0].click()
        sleep(1)
        if pdbid != '':
            # 选择“Single network”选项后，在PDB ID输入框中输入“pdb id”, for example "4BK1"
            self.driver.find_element_by_id('pdbid').send_keys(pdbid.upper())
        if pdb_file != '':
            # 在choose file按钮下上传“pdb file”
            self.driver.find_element_by_id('upload').send_keys(pdb_file)
        if cut_hydro != '':
            # 上传氢键计算的截断值
            input_cut_hydro = self.driver.find_element_by_name('cut_hydro')
            self.driver.execute_script("arguments[0].value = ''", input_cut_hydro)
            input_cut_hydro.send_keys(cut_hydro)
        sleep(2)
        self.driver.find_element_by_name('captcha_code').send_keys(result_four)
        sleep(1)
        self.driver.find_element_by_id('sub1').click()
        sleep(1)
        #识别判断是否提交正确
        try:
            captcha_error = False
            # submit_success = self.driver.find_element_by_xpath("//*[contains(text(),'http://cpclab.uni-duesseldorf.de/cna/results/')]")
            # href = submit_success.get_attribute('textContent')
            urls = self.driver.find_elements_by_xpath("//a")
            for url in urls:
                url_link = url.get_attribute("outerHTML") ##get all the urls
                # url_link = url.get_attribute("href").toString()
                # print(url_link)
                if '/cna/results' in url_link:
                    submit_success = url_link.split(' ')[1].split('=')[-1].strip('"')
                    href = True
                    sleep(1)
                    print(submit_success)
                    self.driver.quit()
                    # url.click()
                    self.get_web_data(submit_success)
                    break
            sleep(1)
            print(href, submit_success)
        except:
            href = False
            submit_success = 'something wrong when submit null'
            # captcha_error = self.driver.find_element_by_xpath("//*[contains(text(),'ERROR: The code you entered was incorrect. Please try again.')]")
            # captcha_error.get_attribute(self)
            sleep(1)
            print(captcha_error)
            self.driver.quit()
        return href, captcha_error, submit_success

if __name__ == '__main__':
##parameters
    work_dir = './'
    os.chdir(work_dir)

    # analysis_type的类型为'm1'--Single network - Single Structure,
    # 'm2'--Ensemble of networks - Ensemble of structures,
    # 'm3'--Ensemble of networks - Single Structure
    analysis_type = 'm1'
    pdbid = ''  ##'5R33'
    # pdb_file = '/home/zhonghua/jinlong/server/stability_analysis/CNA_submit/findholes_workflow_scripts/test_model/5lnw_chainA_schro/1nf2_chainA_delMg.pdb' #'/home/zhonghua/Downloads/PercInd.pdb'
    pdb_file = sys.argv[1]
    cut_hydro = '0.5'  #	Hydrophobic cutoff: Constant value

    if analysis_type == 'm2':
        pdbid = ''  ##必须为输入文件，不能为PDB ID

##main function
    while True:
        while True:
            a = VerificationCode()
            result_four = a.image_str()
            if len(result_four) == 4:
                break
        href, captcha_error, submit_success = a.fill_data(result_four, analysis_type, pdbid, pdb_file, cut_hydro)
        if href is True:
            break
        if captcha_error is True:
            continue
