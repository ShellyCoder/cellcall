# miniconda琐碎记录

### 1. 安装conda
```
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# https://repo.anaconda.com/miniconda/ 所有miniconda的下载地址
# 这个版本是适合于linux的，要看清楚噢。
# 这里选择的是latest-Linux版本，所以下载的程序会随着python的版本更新而更新
#（现在下载的版本默认的python版本已经是3.7了）

chmod 777 Miniconda3-latest-Linux-x86_64.sh #给执行权限
bash Miniconda3-latest-Linux-x86_64.sh #运行

source .bashrc
```


**参考：**[https://www.jianshu.com/p/edaa744ea47d](https://www.jianshu.com/p/edaa744ea47d)
**
### 2. 镜像
#### 2.1 概览
**一般新装的miniconda是没有 .condarc 镜像配置文件的，需要 conda config --set show_channel_urls yes后才会自动生成。如下图可见**
```
conda config --set show_channel_urls yes
```
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1594985099777-f2d8c5a1-8cca-46a4-ad77-7d2b5cc63b87.png#align=left&display=inline&height=287&margin=%5Bobject%20Object%5D&name=image.png&originHeight=574&originWidth=1383&size=115107&status=done&style=none&width=691.5)
** .condarc 镜像配置文件一开始都是默认的配置，因此要加上离我们近的清华，或北师大镜像**
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1594985271784-3a32d89c-66d4-440a-9c61-cf25b8db5482.png#align=left&display=inline&height=181&margin=%5Bobject%20Object%5D&name=image.png&originHeight=362&originWidth=962&size=65427&status=done&style=none&width=481)
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1594985279025-d989ae14-694d-461a-8c5e-1101561d023b.png#align=left&display=inline&height=306&margin=%5Bobject%20Object%5D&name=image.png&originHeight=612&originWidth=1142&size=90176&status=done&style=none&width=571)
即可，配置成功。


#### 2.2 detail

官方channel: (先不要急着添加这两个哦~,只要添加下面的清华的4个镜像地址就足够了的~)：
```
conda config --add channels bioconda
conda config --add channels conda-forge
```

官方的话这两个channel应该就够了的。
`2020-06-14 update：但是其实现在用国内的镜像比较多，官方的频道相较而言速度较慢。但也不是绝对的，有小伙伴跟我说他使用官方的频道也很流畅，所以见仁见智啦。另外，不建议加入大量的相同的频道，如添加了官方的bioconda之后又添加清华的bioconda镜像，没有必要，而且会拖慢速度。`


**清华镜像**
```
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
```
`
2020-06-14 update: 为了分担清华源镜像的压力，最近北京外国语大学也开启了镜像站点，同样是由清华TUNA团队维护的，如果有小伙伴遇到清华源速度很慢的情况的话，可以考虑换成北外的镜像。`
`新闻传送门：[https://mirrors.tuna.tsinghua.edu.cn/news/bfsu-mirror/](https://mirrors.tuna.tsinghua.edu.cn/news/bfsu-mirror/)`
`镜像传送门：[https://mirrors.bfsu.edu.cn/help/anaconda/](https://mirrors.bfsu.edu.cn/help/anaconda/)`
`2020-08-05 update: 为了方便大家(当然主要是自己偷懒用), 把北外的链接也给写出来, 这样就可以直接复制粘贴了~当然两者取其一就可以了, 不用重复添加.

2020-08-10 update: 在生信技能树的群里由群友@合肥-生信-gzcdo 提供了两个新的conda的国内镜像源https://mirrors.nju.edu.cn/anaconda/https://mirrors.sjtug.sjtu.edu.cn/anaconda/各位朋友也可以试试看这两个镜像呀!~`


**北师大镜像**
> channels:
>   - [https://mirrors.bfsu.edu.cn/anaconda/pkgs/main](https://mirrors.bfsu.edu.cn/anaconda/pkgs/main)
>   - [https://mirrors.bfsu.edu.cn/anaconda/pkgs/free](https://mirrors.bfsu.edu.cn/anaconda/pkgs/free)
>   - [https://mirrors.bfsu.edu.cn/anaconda/pkgs/r](https://mirrors.bfsu.edu.cn/anaconda/pkgs/r)
>   - [https://mirrors.bfsu.edu.cn/anaconda/pkgs/pro](https://mirrors.bfsu.edu.cn/anaconda/pkgs/pro)
>   - [https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2](https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2)
>   - [https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge](https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge)
>   - [https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda](https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda)
>   - bioconda
>   - r
>   - conda-forge
>   - defaults


显示安装的频道
```
conda config --set show_channel_urls yes
```
查看已经添加的channels
```
conda config --get channels
```
已添加的channel在哪里查看
```
vim ~/.condarc
```
不希望以开终端base环境默认初始化：
conda config --set auto_activate_base false
### 
3. 创建个人的虚拟环境

conda create -n lty python=3.6
设置好镜像之后速度特别快，1分钟就可以完成。
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1594985446697-00640284-3e9d-4a24-a3eb-78f7c37aece6.png#align=left&display=inline&height=256&margin=%5Bobject%20Object%5D&name=image.png&originHeight=512&originWidth=1460&size=55277&status=done&style=none&width=730)


**查看有哪些虚拟环境**
```
conda env list
```
**![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1594985502717-1c903c5d-af79-4a92-8ff4-611e8b77f8b2.png#align=left&display=inline&height=79&margin=%5Bobject%20Object%5D&name=image.png&originHeight=158&originWidth=1031&size=14549&status=done&style=none&width=515.5)**
**也就我刚创建的环境，还有基础的base环境

**删除也很容易的
```
conda remove -n myenv --all
```
**
### 4. 软件安装和卸载
当然, 也可以用这个命令进行搜索（会稍微慢一点）
```
conda search gatk
```

安装完成后，可以用“which 软件名”来查看该软件安装的位置：
```
which gatk
```

如需要安装特定的版本:
```
conda install 软件名=版本号
conda install gatk=3.7
```

这时conda会先卸载已安装版本，然后重新安装指定版本。
查看已安装软件:
```
conda list
```

更新指定软件:
```
conda update gatk
```

卸载指定软件:
```
conda remove gatk
```


### 5. 创建软件的软链接
跟着命令一路敲到这里的小旁友们估计发现了，现在退出conda环境之后之前安装的软件全都GG了，敲命令没法执行了！
怎么办呢！其实只要把安装好的软件软连接到一个处在环境变量里的位置就可以使用了。三步走：

- 第一步，创建一个文件夹
我一般的习惯是在`/home`目录下创建一个`.soft`文件夹
- 第二步，将这个文件夹添加到环境变量中
```
export PATH="~/.soft:$PATH"
```

- 第三步，软链接
```
ln -s ~/miniconda3/bin/gatk ~/.soft
```
这样就可以运行啦~如果还是不行建议试试初始化一下bashrc：`. ./bashrc`**
**
