## 程序版本管理
## 2017-6-1
## 李扬

.PHONY : veryclean, push, backup

veryclean :
	$(MAKE) -C pwf veryclean
	$(MAKE) -C mtx veryclean
	$(MAKE) -C mtx/test veryclean

backup :
	git add .
	git commit -m "`date` 备份提交"

push : veryclean backup
	git push
