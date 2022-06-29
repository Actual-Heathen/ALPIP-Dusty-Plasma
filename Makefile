TARGET = smallScaleTest

all:
	$(MAKE) -C source
	$(MAKE) clean
clean:
	$(RM) $(TARGET) ./source/$(TARGET)
