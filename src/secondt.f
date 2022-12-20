      function second()
      real*8 second
      real*4 t(2)
c
c     IBM does not have etime(t). Use mclock() instead
c
      sec = mclock()
      second = sec/100.0d0
      return
      end
