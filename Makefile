TARGET = smallScaleTest

all:
	$(MAKE) -C source
clean:
	$(RM) $(TARGET) ./source/$(TARGET)
