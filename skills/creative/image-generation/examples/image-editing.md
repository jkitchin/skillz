# Image Editing Examples

> **Note**: Examples use the `google-genai` SDK. The `google-generativeai` package has been deprecated.

## Overview

Image editing involves modifying existing images using text instructions. This capability is available with Gemini models (both Flash and Pro) and DALL-E 2 (inpainting). This guide shows practical examples of various editing operations.

## Key Capabilities

- **Background replacement** - Change the scene behind a subject
- **Style transfer** - Apply different artistic styles to existing images
- **Object addition** - Add new elements to scenes
- **Color/lighting adjustments** - Modify mood and atmosphere
- **Composition changes** - Alter layout and framing
- **Inpainting** (DALL-E 2) - Replace specific masked regions

## Example 1: Background Replacement

### User Request
"Change the background of this portrait to a sunset beach scene"

### Analysis
- **Operation:** Background replacement
- **Model choice:** Gemini (Flash or Pro depending on quality needs)
- **Input:** Portrait photo with plain background
- **Output:** Same subject with beach sunset background

### Prompt Strategy
Be specific about what stays and what changes.

```
Change the background to a beautiful beach at sunset, golden hour lighting with
orange and pink sky, ocean waves and palm trees visible, keep the subject in
the foreground unchanged, seamless composition, photorealistic blend
```

### Implementation (Gemini)
```python
from google import genai
from google.genai import types
from PIL import Image

client = genai.Client()  # Uses GEMINI_API_KEY env var

# Load existing image
original_image = Image.open("portrait.jpg")

# Edit instruction
instruction = """Change the background to a beautiful beach at sunset, golden
hour lighting with orange and pink sky, ocean waves and palm trees visible in
background, keep the subject in foreground unchanged, seamless composition,
photorealistic blend"""

response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=[original_image, instruction],
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"]
    )
)

# Save edited image
for part in response.parts:
    if part.inline_data is not None:
        image = part.as_image()
        image.save("portrait_beach_background.png")
        print("Saved: portrait_beach_background.png")
```

## Example 2: Style Transfer

### User Request
"Make this photo look like a Van Gogh painting"

### Analysis
- **Operation:** Artistic style transfer
- **Model choice:** Gemini or DALL-E 3
- **Input:** Regular photograph
- **Output:** Same scene in Van Gogh's style

### Prompt
```
Transform this image into the style of Van Gogh's paintings, with characteristic
swirling brushstrokes, vibrant expressive colors, thick impasto texture, post-
impressionist technique, maintain the scene composition but reinterpret in
Van Gogh's artistic style
```

### Implementation (Gemini)
```python
original_image = Image.open("photo.jpg")

instruction = """Transform this into Van Gogh's painting style with swirling
brushstrokes, vibrant expressive colors, thick impasto texture, post-
impressionist technique, maintain composition"""

response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=[original_image, instruction],
    config=types.GenerateContentConfig(response_modalities=["TEXT", "IMAGE"])
)

for part in response.parts:
    if part.inline_data is not None:
        image = part.as_image()
        image.save("photo_vangogh_style.png")
```

## Example 3: Adding Objects

### User Request
"Add a laptop on the desk in this office photo"

### Analysis
- **Operation:** Object addition
- **Model choice:** Gemini models
- **Input:** Office photo without laptop
- **Output:** Same scene with realistic laptop added

### Prompt
```
Add a modern laptop computer on the desk, silver MacBook-style design, screen
facing viewer showing a generic desktop, positioned naturally on the desk
surface, realistic lighting and shadows matching the scene, seamless integration
```

### Implementation
```python
office_photo = Image.open("office_desk.jpg")

instruction = """Add a modern laptop computer on the desk, silver MacBook-style,
screen showing generic desktop, positioned naturally with realistic lighting
and shadows, seamless integration"""

response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=[office_photo, instruction],
    config=types.GenerateContentConfig(response_modalities=["TEXT", "IMAGE"])
)
```

## Example 4: Lighting and Mood Changes

### User Request
"Make this daytime scene look like night with city lights"

### Analysis
- **Operation:** Lighting transformation
- **Model choice:** Gemini Pro (complex lighting changes)
- **Input:** Daytime cityscape
- **Output:** Night version with artificial lighting

### Prompt
```
Transform this to a night scene, dark blue evening sky with stars visible,
building windows glowing with warm interior lights, street lamps illuminated,
car headlights and taillights, neon signs lit up, atmospheric night mood,
maintain architecture and composition
```

### Implementation (Gemini Pro)
```python
cityscape = Image.open("city_day.jpg")

instruction = """Transform to night scene, dark blue sky with stars, building
windows glowing with warm light, street lamps on, car lights, neon signs,
atmospheric night mood, keep architecture intact"""

response = client.models.generate_content(
    model="gemini-3-pro-image-preview",
    contents=[cityscape, instruction],
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        image_config=types.ImageConfig(image_size="2K")
    )
)
```

## Example 5: Color Grading

### User Request
"Give this photo a warm vintage film look"

### Analysis
- **Operation:** Color grading and filter application
- **Model choice:** Gemini Flash (adequate for color changes)
- **Input:** Modern digital photo
- **Output:** Vintage-styled version

### Prompt
```
Apply warm vintage film color grading, slightly faded colors with warm yellow-
orange tones, subtle grain texture, slightly reduced contrast, 1970s film
photography aesthetic, nostalgic mood, maintain subject and composition
```

### Implementation
```python
photo = Image.open("modern_photo.jpg")

instruction = """Apply warm vintage film color grading, faded colors with warm
yellow-orange tones, subtle grain, reduced contrast, 1970s film aesthetic,
nostalgic mood"""

response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=[photo, instruction],
    config=types.GenerateContentConfig(response_modalities=["TEXT", "IMAGE"])
)
```

## Example 6: Season Change

### User Request
"Make this summer landscape look like winter with snow"

### Analysis
- **Operation:** Seasonal transformation
- **Model choice:** Gemini Pro (complex environmental changes)
- **Input:** Summer landscape
- **Output:** Winter version with snow

### Prompt
```
Transform this to a winter scene, add snow covering the ground and trees,
bare tree branches with snow, overcast winter sky, cold blue color palette,
frost on surfaces, winter atmosphere, maintain landscape composition and
major features
```

### Implementation
```python
summer_landscape = Image.open("summer_scene.jpg")

instruction = """Transform to winter scene with snow covering ground and trees,
bare snow-covered branches, overcast winter sky, cold blue tones, frost,
winter atmosphere, keep landscape composition"""

response = client.models.generate_content(
    model="gemini-3-pro-image-preview",
    contents=[summer_landscape, instruction],
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        image_config=types.ImageConfig(image_size="2K")
    )
)
```

## Example 7: DALL-E 2 Inpainting

### User Request
"Replace the person's shirt with a different color and style"

### Analysis
- **Operation:** Specific region replacement
- **Model choice:** DALL-E 2 (supports masked editing)
- **Input:** Photo + mask indicating shirt area
- **Output:** Same photo with modified shirt

### Process
1. Create a mask (transparent PNG where transparent = area to change)
2. Provide mask and instruction to DALL-E 2
3. Generate edited version

### Implementation (DALL-E 2)
```python
from openai import OpenAI
import requests
from pathlib import Path

client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

# Original image and mask (mask is transparent PNG with shirt area transparent)
with open("person_photo.png", "rb") as image_file, \
     open("shirt_mask.png", "rb") as mask_file:

    instruction = """Professional blue button-down shirt on the person,
    well-fitted, professional appearance, realistic fabric texture and lighting"""

    response = client.images.edit(
        model="dall-e-2",
        image=image_file,
        mask=mask_file,
        prompt=instruction,
        n=1,
        size="1024x1024"
    )

    edited_url = response.data[0].url
    edited_image = requests.get(edited_url).content
    Path("person_new_shirt.png").write_bytes(edited_image)
```

### Creating Masks
```python
# Example: Create mask programmatically with PIL
from PIL import Image, ImageDraw

# Load original
img = Image.open("person_photo.png").convert("RGBA")

# Create mask (white background, transparent where you want to edit)
mask = Image.new("RGBA", img.size, (255, 255, 255, 255))
draw = ImageDraw.Draw(mask)

# Define the region to edit (shirt area) - make it transparent
# This example uses a simple rectangle, but you'd want precise selection
draw.ellipse([200, 300, 600, 800], fill=(0, 0, 0, 0))  # Transparent

mask.save("shirt_mask.png")
```

## Example 8: Remove and Replace Objects

### User Request
"Remove the car from this street scene and replace it with a bicycle"

### Analysis
- **Operation:** Object removal and replacement
- **Model choice:** Gemini models
- **Input:** Street with car
- **Output:** Same street with bicycle instead

### Prompt
```
Remove the car from the scene and replace it with a vintage bicycle leaning
against the lamppost, maintain the street scene composition, realistic
integration with appropriate scale and lighting, seamless edit
```

### Implementation
```python
street_scene = Image.open("street_with_car.jpg")

instruction = """Remove the car and replace with vintage bicycle leaning
against lamppost, maintain composition, realistic scale and lighting,
seamless edit"""

response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=[street_scene, instruction],
    config=types.GenerateContentConfig(response_modalities=["TEXT", "IMAGE"])
)
```

## Example 9: Perspective and Composition Changes

### User Request
"Change the angle to make it look like it was shot from above"

### Analysis
- **Operation:** Perspective transformation
- **Model choice:** Gemini Pro (complex spatial understanding)
- **Input:** Eye-level photo
- **Output:** Overhead/bird's eye perspective

### Prompt
```
Transform this to an overhead bird's eye view perspective, looking down from
above, maintain the subject and scene elements but recompose from aerial
viewpoint, realistic spatial relationships, professional overhead photography
```

### Implementation
```python
eye_level_photo = Image.open("scene_eye_level.jpg")

instruction = """Transform to overhead bird's eye view, looking down from above,
maintain elements but recompose from aerial viewpoint, realistic spatial
relationships"""

response = client.models.generate_content(
    model="gemini-3-pro-image-preview",
    contents=[eye_level_photo, instruction],
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        image_config=types.ImageConfig(image_size="2K")
    )
)
```

## Example 10: Facial Expression Changes

### User Request
"Make the person in this photo smiling"

### Analysis
- **Operation:** Facial feature modification
- **Model choice:** Gemini models
- **Input:** Portrait with neutral expression
- **Output:** Same person smiling

### Prompt
```
Change the person's facial expression to a genuine warm smile, natural and
friendly expression, keep all other features and details identical, maintain
photorealistic appearance, seamless edit
```

### Implementation
```python
portrait = Image.open("neutral_expression.jpg")

instruction = """Change facial expression to genuine warm smile, natural and
friendly, keep all other features identical, photorealistic, seamless"""

response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=[portrait, instruction],
    config=types.GenerateContentConfig(response_modalities=["TEXT", "IMAGE"])
)
```

## Multi-Step Iterative Editing

### User Request
"First change the background, then add a sunset effect, then make it warmer"

### Analysis
- **Operation:** Sequential edits building on each other
- **Model choice:** Gemini with chat interface
- **Approach:** Use conversation context for iterations

### Implementation (Gemini Chat)
```python
# Create chat for multi-turn editing
chat = client.chats.create(
    model="gemini-2.5-flash-image",
    config=types.GenerateContentConfig(response_modalities=["TEXT", "IMAGE"])
)

# Load original
original = Image.open("portrait.jpg")

# Edit 1: Background change
response1 = chat.send_message([original, "Change the background to a city skyline at dusk"])

for part in response1.parts:
    if part.inline_data is not None:
        image = part.as_image()
        image.save("edit_step1.png")

# Edit 2: Add sunset lighting (builds on previous)
response2 = chat.send_message("Add dramatic sunset lighting with golden rays coming from the left side")

for part in response2.parts:
    if part.inline_data is not None:
        image = part.as_image()
        image.save("edit_step2.png")

# Edit 3: Warm color grading
response3 = chat.send_message("Apply warm color grading with enhanced orange and yellow tones")

for part in response3.parts:
    if part.inline_data is not None:
        image = part.as_image()
        image.save("edit_final.png")
```

## Best Practices for Image Editing

### 1. Be Specific About What Changes
```
Bad: "Make it better"
Good: "Increase the brightness by 20%, add more contrast, enhance the blue tones in the sky"
```

### 2. Specify What Should Stay the Same
```
"Change the background to a forest scene, but keep the subject in the foreground
exactly as they are, including clothing, pose, and lighting on the subject"
```

### 3. Describe Integration
```
"Add realistic shadows and lighting that match the existing scene, seamless
integration, maintain consistent perspective and scale"
```

### 4. Use Reference Terms
```
"Apply the same lighting style as the original", "Match the existing color
palette", "Maintain the current composition and framing"
```

### 5. Layer Complex Edits
For major changes, break into steps:
1. First major structural change
2. Then lighting/color adjustments
3. Finally detail refinements

### 6. Choose Right Model
- **Simple edits** (color, brightness) → Gemini Flash
- **Complex transformations** (perspective, major additions) → Gemini Pro
- **Specific masked regions** → DALL-E 2 inpainting
- **Quick iterations** → Start with Flash, upgrade if needed

### 7. Maintain Context with Chat
For multi-step edits, use Gemini's chat interface to maintain context and build progressively.

## Common Editing Patterns

### Background Replacement
```
[original_image] + "Change the background to [new background description],
keep the subject unchanged, seamless composition, photorealistic blend"
```

### Style Transfer
```
[original_image] + "Apply [art style] to this image, maintain composition
and subject but reinterpret in [style] technique"
```

### Object Addition
```
[original_image] + "Add [object description] to [location in scene],
realistic integration with matching lighting and shadows"
```

### Lighting Change
```
[original_image] + "Transform the lighting to [lighting description],
maintain subject and composition but change atmosphere to [mood]"
```

### Color Grading
```
[original_image] + "Apply [color grade description], [tone description],
maintain subject details and composition"
```

## Troubleshooting

**Edit doesn't look realistic:**
- Add "seamless integration", "realistic lighting", "natural blend" to prompt
- Use higher quality model (Gemini Pro)
- Try more specific instructions

**Subject changed when only background should:**
- Explicitly state "keep subject identical/unchanged"
- Use more specific description of what should remain

**Colors look off:**
- Specify "match existing lighting conditions"
- Add "maintain color harmony" or "consistent color palette"

**Multiple elements changed unexpectedly:**
- Be more specific about what should and shouldn't change
- Break into multiple steps with chat interface

## Summary

Image editing with AI is powerful for:
- Background replacement
- Style transfers
- Adding/removing objects
- Lighting and mood changes
- Color grading and filters
- Perspective transformations

**Key to success:** Be specific about both what changes and what stays the same, describe integration requirements, and use the appropriate model for the complexity of the edit.
